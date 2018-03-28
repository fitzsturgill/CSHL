% multilinear regression with licking


[filename, filepath, ind] = uigetfile('TE.mat', 'select TE to load', 'MultiSelect', 'off');


%%
channel = 2;

load(fullfile(filepath, filename));
disp(['*** Loaded ' fullfile(filepath, filename)]);

savepath = fullfile(filepath, ['lickRegression_ch' num2str(channel)]);
ensureDirectory(savepath);
saveOn = 1;
%%
TE.licksTS = bpEventToTimeSeries(TE, 'Port1In', 'duration', 11);

%
% Initial Conditions
Fs = 20; % sample rate

Kw = [-0.5 1] ; % Kernel/weight offset window
% Kw = [0 2] ; % Kernel/weight offset window


% generate weight time offset indices
Kw2 = Kw * Fs;
Kwix = Kw2(1):Kw2(2);
margins = [sum(Kwix < 0) sum(Kwix > 0)]; % to avoid edge effects due to using circ-shifted predictor variables
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
DMraw = zeros(truncSize(1), truncSize(2), length(Kwix));


% shift matrices
for counter = 1:length(Kwix)
    shifted = circshift(dm, Kwix(counter), 1);
    DMraw(:,:,counter) = shifted(margins(1) + 1 : end - margins(2), :);
end

% shuffle and/or reshape
DMshuff = DMraw(:,randperm(size(DMraw, 2)), :);
DM = reshape(DMraw, prod(truncSize), length(Kwix));
DMshuff = reshape(DMshuff, prod(truncSize), length(Kwix));
% try adding on time
% DM = [DM repmat(diag(ones(truncSize(1), 1)), truncSize(2), 1)]; 
% DMshuff = [DMshuff repmat(diag(ones(truncSize(1), 1)), truncSize(2), 1)]; 
% regress
[w wconf R] = regress(T, DM);
% cross validate
Rsq = regressCrossValidate(T, DM);

nShuffles = 50;
wShuff_all = zeros(kSize, nShuffles);
R_all = zeros(length(R), nShuffles);

for counter = 1:nShuffles
    disp(['Shuffle # ' num2str(counter)]);
    DMshuff = DMraw(:,randperm(size(DMraw, 2)), :);
    DMshuff = reshape(DMshuff, prod(truncSize), length(Kwix));
    [wshuff wshuff_sem R_shuff] = regress(T, DMshuff);
%     Rsq_shuff = regressCrossValidate(T, DMshuff);
    wShuff_all(:,counter) = wshuff;
    R_all(:,counter) = R_shuff;
end

%
R_shuff = mean(R_all, 2);
wshuff = mean(wShuff_all, 2);
wshuff_sem = std(wShuff_all, 0, 2) / sqrt(size(wShuff_all, 2));
w_corrected = w - wshuff;

prediction_corrected = DM * w_corrected;
R_corrected = T - prediction_corrected;


%
% xData for plotting weight kernel
% xData_time = linspace(0, 11, truncSize(1));
xData = linspace(Kw(1), Kw(2), kSize);
%
% ensureFigure('test', 1); boundedline(xData', w, wconf - w);
% ensureFigure('test', 1); plot(xData', w, 'b'); hold on; plot(xData', wconf,'c')
R = reshape(R, truncSize(1), truncSize(2));
R = R';
T = reshape(T, truncSize(1), truncSize(2));
T = T';
R_corrected = reshape(R_corrected, truncSize(1), truncSize(2));
R_corrected = R_corrected';
R_shuff = reshape(R_shuff, truncSize(1), truncSize(2));
R_shuff = R_shuff';

%%
savename = ['LickModel_ch' num2str(channel)];
ensureFigure(savename, 1);
subplot(2,2,1); imagesc(T, [-2 2]); title('target');
subplot(2,2,2); imagesc(R_corrected, [-2 2]); title('residuals');
subplot(2,2,3); 
% plot(xData', w(1:kSize), 'b'); hold on; plot(xData', wconf(1:kSize, :),'c'); 
boundedline(xData', w, wconf(:,1) - w, 'b'); hold on;
xlabel('time (s)'); ylabel('Beta weights'); title('Kernel');
boundedline(xData', wshuff, wshuff_sem, 'k'); %plot(xData', wconf_shuff(1:kSize, :), 'm');
subplot(2,2,4); plot(xData', w - wshuff); title('Kernel corrected');

saveas(gcf, fullfile(savepath, savename), 'fig');
saveas(gcf, fullfile(savepath, savename), 'jpeg');
saveas(gcf, fullfile(savepath, savename), 'epsc');

%% save
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_39\lick_regression\';

licks = struct(...
    'weights_licks', w(1:kSize),...
    'wconf_licks', wconf(1:kSize,:),...    
    'R', R,...
    'Rsq', Rsq,...
    'xData', xData,...
    'includedTrials', include...
    );

licks_shuff = struct(...
    'weights_licks', wshuff(1:kSize),...
    'wconf_licks', wshuff_sem(1:kSize,:),...
    'R', R_shuff,...
    'Rsq', Rsq_shuff,...
    'xData', xData,...
    'includedTrials', include...
    );
licks_corrected = struct(...
    'weights_licks', w_corrected(1:kSize),...
    'R', R_corrected,...
    'xData', xData,...
    'includedTrials', include...
    );

save(fullfile(savepath, 'licks.mat'), 'licks*');


%% cuedOutcome plot

savename = 'cuedOutcome_regress_averages';
ensureFigure(savename, 1);
cuedOutcome_Conditions;

subplot(2,1,1); 
rXData = TE.Photometry.xData(margins(1)+1:end - margins(2));
plot(rXData, mean(licks_corrected.R(highValueTrials & include & rewardTrials, :)), 'g'); hold on;
plot(rXData, mean(licks_corrected.R(lowValueTrials & include & rewardTrials, :)), 'm');
plot(rXData, mean(licks_corrected.R(uncuedTrials & include & rewardTrials, :)), 'k');
set(gca, 'XLim', [rXData(1) rXData(end)]);
title('Residuals');

subplot(2,1,2);
[ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1, 'FluorDataField', 'ZS', 'linespec', {'g', 'm', 'k'}); %high value, reward  
title('reward conditions'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 
set(gca, 'XLim', [rXData(1) rXData(end)]);
legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');

saveas(gcf, fullfile(savepath, savename), 'fig');
saveas(gcf, fullfile(savepath, savename), 'jpeg');
saveas(gcf, fullfile(savepath, savename), 'epsc');

%% LNL Odor plot

savename = ['LNL_regress_averages_ch' num2str(channel)'];
ensureFigure(savename, 1);
LNL_conditions;

subplot(2,1,1); 
rXData = TE.Photometry.xData(margins(1)+1:end - margins(2));
plot(rXData, mean(licks_corrected.R(csPlusTrials & hitTrials & include & rewardTrials, :)), 'c'); hold on;
plot(rXData, mean(licks_corrected.R(csPlusTrials & missTrials & include & rewardTrials, :)), 'm');
plot(rXData, mean(licks_corrected.R(include & uncuedTrials & rewardTrials, :)), 'k');
set(gca, 'XLim', [rXData(1) rXData(end)]);
title('Residuals');

trialTypesLNL = {csPlusTrials & hitTrials & include & rewardTrials, csPlusTrials & missTrials & include & rewardTrials,...
    include & uncuedTrials & rewardTrials};
subplot(2,1,2);
[ha, hl] = phPlotAverageFromTE(TE, trialTypesLNL, channel, 'FluorDataField', 'ZS', 'linespec', {'c', 'm', 'k'}); %high value, reward  
title('reward conditions'); ylabel('Z Score'); xlabel('time from cue(s)'); 
set(gca, 'XLim', [rXData(1) rXData(end)]);
legend(hl, {'CS+, hit', 'CS+, miss', 'uncued, rew'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
% 
% saveas(gcf, fullfile(savepath, savename), 'fig');
% saveas(gcf, fullfile(savepath, savename), 'jpeg');
% saveas(gcf, fullfile(savepath, savename), 'epsc');



