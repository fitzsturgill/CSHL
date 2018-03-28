% multilinear regression with licking

% load('Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_37\TE.mat');
% load('Z:\SummaryAnalyses\LickNoLick_odor_v2\DC_35\TE.mat');
% load('Z:\SummaryAnalyses\LickNoLick_odor_v2\DC_36\TE.mat'); % best
% load('Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_36\TE.mat');
% load('Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_35\TE.mat');
% load('Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_25\TE.mat');
load('Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_39\TE.mat');
% LNL_conditions;
%
TE.licksTS = bpEventToTimeSeries(TE, 'Port1In', 'duration', 11);

%
% Initial Conditions
Fs = 20; % sample rate
channel = 1;
Kw = [-1 1]; % Kernel/weight offset window


% generate weight time offset indices
Kw2 = Kw * Fs;
Kwix = Kw2(1):Kw2(2);
margins = [sum(Kwix < 0) sum(Kwix) > 0]; % to avoid edge effects due to using circ-shifted predictor variables
kSize = length(Kwix);

% Target
% find any NaN containing trials and exclude them
temp = TE.Photometry.data(channel).ZS;
include = sum(isnan(temp), 2) == 0; 
T = temp(include, :);
T=T';
T = T(margins(1) + 1 : end - margins(2), :);
truncSize = size(T);
% T = ([1 1 1 1; 2 2 2 2; 3 3 3 3]);
T = T(:);
T = zscore(T);

% Design matrix/ Predictors
dm = TE.licksTS.data(include, :)' > 0; % get rid of spurious fast bursts of licks
% dm = ([0:100; 50:150; 100:200])';
DM = zeros(truncSize(1), truncSize(2), length(Kwix));


% shift matrices
for counter = 1:length(Kwix)
    shifted = circshift(dm, Kwix(counter), 1);
    DM(:,:,counter) = shifted(margins(1) + 1 : end - margins(2), :);
end

% shuffle and/or reshape
DMshuff = DM(:,randperm(size(DM, 2)), :);
DM = reshape(DM, prod(truncSize), length(Kwix));
DMshuff = reshape(DMshuff, prod(truncSize), length(Kwix));
% try adding on time
DM = [DM repmat(diag(ones(truncSize(1), 1)), truncSize(2), 1)]; 
DMshuff = [DMshuff repmat(diag(ones(truncSize(1), 1)), truncSize(2), 1)]; 
% regress
[w wconf R] = regress(T, DM);
% cross validate
Rsq = regressCrossValidate(T, DM);

[wshuff wconf_shuff R_shuff] = regress(T, DMshuff);
Rsq_shuff = regressCrossValidate(T, DMshuff);
%
% xData for plotting weight kernel
xData_time = linspace(0, 11, truncSize(1));
xData = linspace(Kw(1), Kw(2), kSize);
%
% ensureFigure('test', 1); boundedline(xData', w, wconf - w);
% ensureFigure('test', 1); plot(xData', w, 'b'); hold on; plot(xData', wconf,'c')
R = reshape(R, truncSize(1), truncSize(2));
R = R';
T = reshape(T, truncSize(1), truncSize(2));
T = T';
%%
ensureFigure(['Lick_Regression' num2str(Kw(1)) 'to' num2str(Kw(2))], 1); subplot(2,2,1); imagesc(T, [-2 2]); title('target');
subplot(2,2,2); imagesc(R, [-2 2]); title('residuals');
subplot(2,2,3); plot(xData', w(1:kSize), 'b'); hold on; plot(xData', wconf(1:kSize, :),'c'); xlabel('time (s)'); ylabel('Beta weights'); title('Kernel');
plot(xData', wshuff(1:kSize), 'r'); plot(xData', wconf_shuff(1:kSize, :), 'm');
 subplot(2,2,4); plot(xData_time, w(kSize + 1:end));
% subplot(2,2,4); plot(xData', w - wshuff);

%% save
savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_39\lick_regression\';

time_and_licks = struct(...
    'weights_licks', w(1:kSize),...
    'wconf_licks', wconf(1:kSize,:),...
    'weights_time', w(kSize+1:end),...
    'wconf_time', wconf(kSize+1:end,:),...    
    'R', R,...
    'Rsq', Rsq,...
    'xData', xData,...
    'xData_time', xData_time,...
    'includedTrials', include...
    );

time_and_licks_shuff = struct(...
    'weights_licks', wshuff(1:kSize),...
    'wconf_licks', wconf_shuff(1:kSize,:),...
    'weights_time', wshuff(kSize+1:end),...
    'wconf_time', wconf_shuff(kSize+1:end,:),...    
    'R', R_shuff,...
    'Rsq', Rsq_shuff,...
    'xData', xData,...
    'xData_time', xData_time,...
    'includedTrials', include...
    );

save(fullfile(savepath, 'time_and_licks.mat'), 'time_and_licks*');


%% plot vs ChAT_39 old data

figure; plot(time_and_licks.xData_time, mean(time_and_licks.R(highValueTrials & include & rewardTrials, :)), 'g'); hold on;
plot(time_and_licks.xData_time, mean(time_and_licks.R(lowValueTrials & include & rewardTrials, :)), 'm');



