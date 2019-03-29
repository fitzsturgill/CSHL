saveOn = 0;
%%
saveOn = 1;
%%
sessions = bpLoadSessions; % load sessions
%% 
TE = makeTE_wheel_alternateLEDs(sessions); % make TE structure

%% Now saved in directory according to first session filename
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = 'Z:\SummaryAnalyses\wheel_405test\';
if length(unique(TE.filename)) > 1
    sep = strfind(TE.filename{1}, '_');
    subjectName = TE.filename{1}(1:sep(2)-1);
else
    sep = strfind(TE.filename{1}, '.');
    subjectName = TE.filename{1}(1:sep(1)-1);
end
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);


%%

dFFMode = {'simple', 'simple'};
try 
    bl = [0 sessions(1).SessionData.NidaqData{1,2}.duration(1) - 0.2];
catch
    bl = [0 9.9];
end

% channel 1 demodulated via channel 1 reference
TE.Photometry_ch1 = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', 1, 'baseline', bl, 'startField', 'Baseline', 'downsample', 305, 'forceAmp', 1);

% channel 1 demodulated via channel 2 reference
TE.Photometry_ch2 = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', 1, 'baseline', bl, 'startField', 'Baseline', 'downsample', 305, 'refChannels', 2, 'forceAmp', 1);


if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%%

LED1trials = TE.LED1_amp > 0 & TE.LED2_amp == 0;
LED2trials = TE.LED2_amp > 0 & TE.LED1_amp == 0;
LED12trials = TE.LED1_amp > 0 & TE.LED2_amp > 0;

%% do bleach fit manually on a trial by trial basis
nTrials = length(TE.filename);
nPoints = size(TE.Photometry_ch1.data(1).raw, 2);
TE.Photometry_ch1.data(1).bleach_fit = zeros(nTrials, nPoints);
TE.Photometry_ch2.data(1).bleach_fit = zeros(nTrials, nPoints);
TE.Photometry_ch1.data(1).bleach_dF = zeros(nTrials, nPoints);
TE.Photometry_ch2.data(1).bleach_dF = zeros(nTrials, nPoints);

x = (0:nPoints - 1);
startFitTime = 1; % start bleach fit n seconds into each trial
Fs = TE.Photometry_ch1.sampleRate;
for counter = 1:nTrials
    if TE.LED1_amp(counter)
        trialData = TE.Photometry_ch1.data(1).raw(counter, :);
        if any(isnan(trialData))
            trialData = inpaint_nans(trialData);
        end
        fo = fitoptions('Method', 'NonlinearLeastSquares',...
            'Upper', [Inf Inf 1/bpX2pnt(1,20)],...
            'Lower', [0 -Inf -1/bpX2pnt(1,20)],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
            'StartPoint', [min(trialData) range(trialData)/2 -1/bpX2pnt(10,20)]...
            );
        model = 'a + b*exp(c*x)';
        ft = fittype(model, 'options', fo);        
        [fitobject, gof, output] = ... % new
            fit(x(bpX2pnt(startFitTime, Fs):end)', trialData(bpX2pnt(startFitTime, Fs):end)', ft, fo);        
        trialFit = fitobject.a + fitobject.b * exp(fitobject.c * x);
        TE.Photometry_ch1.data(1).bleach_fit(counter,:) = trialFit;
        TE.Photometry_ch1.data(1).bleach_dF(counter,:) = trialData - trialFit;
    end
    
    if TE.LED2_amp(counter)
        trialData = TE.Photometry_ch2.data(1).raw(counter, :);
        if any(isnan(trialData))
            trialData = inpaint_nans(trialData);
        end        
        fo = fitoptions('Method', 'NonlinearLeastSquares',...
            'Upper', [Inf Inf 1/bpX2pnt(1,20)],...
            'Lower', [0 -Inf -1/bpX2pnt(1,20)],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
            'StartPoint', [min(trialData) range(trialData)/2 -1/bpX2pnt(10,20)]...
            );
        model = 'a + b*exp(c*x)';
        ft = fittype(model, 'options', fo);                
        [fitobject, gof, output] = ... % new
            fit(x(bpX2pnt(startFitTime, Fs):end)', trialData(bpX2pnt(startFitTime, Fs):end)', ft, fo);        
        trialFit = fitobject.a + fitobject.b * exp(fitobject.c * x);
        TE.Photometry_ch2.data(1).bleach_fit(counter,:) = trialFit;
        TE.Photometry_ch2.data(1).bleach_dF(counter,:) = trialData - trialFit;
    end    
end
%% some checks of bleach fits
saveName = 'bleach_check';
ensureFigure(saveName, 1);
subplot(3,4,1); imagesc(TE.Photometry_ch1.data(1).bleach_fit(LED12trials, :)); ylabel('bleach fit'); title('470, both LEDs');
subplot(3,4,5); imagesc(TE.Photometry_ch1.data(1).bleach_dF(LED12trials, :)); ylabel('bleach corrected');
subplot(3,4,9); imagesc(TE.Photometry_ch1.data(1).raw(LED12trials, :)); ylabel('raw');
subplot(3,4,2); imagesc(TE.Photometry_ch1.data(1).bleach_fit(LED1trials, :)); title('470 only');
subplot(3,4,6); imagesc(TE.Photometry_ch1.data(1).bleach_dF(LED1trials, :)); colorbar;
subplot(3,4,10); imagesc(TE.Photometry_ch1.data(1).raw(LED1trials, :)); colorbar;
subplot(3,4,3); imagesc(TE.Photometry_ch2.data(1).bleach_fit(LED12trials, :)); title('405, both LEDs');
subplot(3,4,7); imagesc(TE.Photometry_ch2.data(1).bleach_dF(LED12trials, :));
subplot(3,4,11); imagesc(TE.Photometry_ch2.data(1).raw(LED12trials, :));%, [-1 1]);
subplot(3,4,4); imagesc(TE.Photometry_ch2.data(1).bleach_fit(LED2trials, :)); title('405 only');
subplot(3,4,8); imagesc(TE.Photometry_ch2.data(1).bleach_dF(LED2trials, :));
subplot(3,4,12); imagesc(TE.Photometry_ch2.data(1).raw(LED2trials, :));%, [-1 1]);

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
end

%% do 405 correction 

% optoin A: use pre-reward period
% preWindow = [-2 0];
% [pre_ch1, ts, tn] = extractDataByTimeStamps(TE.Photometry_ch1.data(1).bleach_dF, TE.Photometry_ch1.startTime, 20, TE.Reward, preWindow, LED12trials);
% [pre_ch2, ts, tn] = extractDataByTimeStamps(TE.Photometry_ch2.data(1).bleach_dF, TE.Photometry_ch2.startTime, 20, TE.Reward, preWindow, LED12trials);
% xData = pre_ch2(find(isfinite(pre_ch2)));
% yData = pre_ch1(find(isfinite(pre_ch1)));

% option B: using entire data set minus the first 2 seconds of each acquisition
xData = TE.Photometry_ch2.data(1).bleach_dF(LED12trials,40:end); xData = xData(:);
yData = TE.Photometry_ch1.data(1).bleach_dF(LED12trials,40:end); yData = yData(:);
saveName = 'linearFit';
ensureFigure(saveName, 1);
scatter(xData, yData, 14, [0 0 1], '.'); hold on;
fo = fitoptions('poly1', 'Robust', 'Bisquare');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
fob = fit(xData, yData, 'poly1', fo); 
fph=plot(fob, 'predfunc'); legend off;
set(fph, 'LineWidth', 0.5, 'Color', [1 0 0]);
textBox(sprintf('470estimated = %.2g  *  405measured + %.2g', fob.p1, fob.p2));
xlabel('405nm'); ylabel('470nm');
title('linear regression using baseline period');
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
end
TE.Photometry_ch1.data(1).dF_corrected = NaN(size(TE.Photometry_ch1.data(1).raw)); % initialize
TE.Photometry_ch1.data(1).dF_estimated = NaN(size(TE.Photometry_ch1.data(1).raw)); % initialize
TE.Photometry_ch1.data(1).dF_estimated(LED12trials, :) = TE.Photometry_ch2.data(1).bleach_dF(LED12trials, :) .* fob.p1 + fob.p2;
TE.Photometry_ch1.data(1).dF_corrected(LED12trials, :) = TE.Photometry_ch1.data(1).bleach_dF(LED12trials, :) - TE.Photometry_ch1.data(1).dF_estimated(LED12trials, :);

if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%% plot some example trials

saveName = 'examples_correction';
ensureFigure(saveName, 1);
trials = find(LED12trials, 4);
for counter = 1:4
    trial = trials(counter);
    subplot(2,2,counter);
    plot(TE.Photometry_ch1.xData, TE.Photometry_ch1.data.bleach_dF(trial, :)', 'r'); hold on;
    plot(TE.Photometry_ch1.xData, TE.Photometry_ch1.data.dF_estimated(trial, :)', 'b');
    plot(TE.Photometry_ch1.xData, TE.Photometry_ch1.data.dF_corrected(trial, :)', 'g');
    legend({'raw', 'estimated', 'corrected'});
end


if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
end

%%
window = [-4 4];
fs = 20; % sample rate
blSamples = (0 - window(1)) * fs;
fluorField = 'bleach_dF';

% LED 470 and 404
[rewards_ch1_both, ts, tn] = extractDataByTimeStamps(TE.Photometry_ch1.data(1).dF_corrected, TE.Photometry_ch1.startTime, 20, TE.Reward, window, LED12trials);
rewards_ch2_both = extractDataByTimeStamps(TE.Photometry_ch2.data(1).(fluorField), TE.Photometry_ch1.startTime, 20, TE.Reward, window, LED12trials);
bl_ch1_both = nanmean(rewards_ch1_both(:,1:blSamples), 2);
rewards_ch1_both = bsxfun(@minus, rewards_ch1_both, bl_ch1_both);
bl_ch2_both = nanmean(rewards_ch2_both(:,1:blSamples), 2);
rewards_ch2_both = bsxfun(@minus, rewards_ch2_both, bl_ch2_both);

% LED 470 only
[rewards_ch1_470, ts, tn] = extractDataByTimeStamps(TE.Photometry_ch1.data(1).(fluorField), TE.Photometry_ch2.startTime, 20, TE.Reward, window, LED1trials);
rewards_ch2_470 = extractDataByTimeStamps(TE.Photometry_ch2.data(1).raw, TE.Photometry_ch2.startTime, 20, TE.Reward, window, LED1trials);
bl_ch1_470 = nanmean(rewards_ch1_470(:,1:blSamples), 2);
rewards_ch1_470 = bsxfun(@minus, rewards_ch1_470, bl_ch1_470);
bl_ch2_470 = nanmean(rewards_ch2_470(:,1:blSamples), 2);
rewards_ch2_470 = bsxfun(@minus, rewards_ch2_470, bl_ch2_470);

% LED 405 only
[rewards_ch1_405, ts, tn] = extractDataByTimeStamps(TE.Photometry_ch1.data(1).(fluorField), TE.Photometry_ch1.startTime, 20, TE.Reward, window, LED2trials);
rewards_ch2_405 = extractDataByTimeStamps(TE.Photometry_ch2.data(1).(fluorField), TE.Photometry_ch2.startTime, 20, TE.Reward, window, LED2trials);
bl_ch1_405 = nanmean(rewards_ch1_405(:,1:blSamples), 2);
rewards_ch1_405 = bsxfun(@minus, rewards_ch1_405, bl_ch1_405);
bl_ch2_405 = nanmean(rewards_ch2_405(:,1:blSamples), 2);
rewards_ch2_405 = bsxfun(@minus, rewards_ch2_405, bl_ch2_405);


% calc mean and SD for ch1 and ch2 (both LEDs on), use this across all LED
% conditions
sd_ch1 = nanmean(nanstd(rewards_ch1_both(:,1:blSamples)));
sd_ch2 = nanmean(nanstd(rewards_ch2_both(:,1:blSamples)));

%% stitch together 
% nTrials = length(TE.filename);
% nSamples = size(rewards_dat_12, 2);
% rewards_dat_all = NaN(nTrials, nSamples);
% rewards_chat_all = NaN(nTrials, nSamples);
% 
% rewards_dat_all(LED12trials, :) = rewards_dat_12(:,:);
% rewards_dat_all(LED1trials, :) = rewards_dat_1(:,:);
% rewards_dat_all(LED2trials, :) = rewards_dat_2(:,:);
% 
% rewards_chat_all(LED12trials, :) = rewards_chat_12(:,:);
% rewards_chat_all(LED1trials, :) = rewards_chat_1(:,:);
% rewards_chat_all(LED2trials, :) = rewards_chat_2(:,:);



climFactor = 4;
clim_ch1 = [-climFactor * sd_ch1, climFactor * sd_ch1];
clim_ch2 = [-climFactor * sd_ch2, climFactor * sd_ch2];

ensureFigure('random_rewards_alternateLEDs', 1); 

% first column, ch1
subplot(3,2,1); image(rewards_ch1_both, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch1); colormap('jet');  ylabel('deltaF, both LEDs'); title('Ch 1');


subplot(3,2,3); image(rewards_ch1_470, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch1); colormap('jet');  ylabel('LED1 only');


subplot(3,2,5); image(rewards_ch1_405, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch1); colormap('jet');  ylabel('LED2 only');

% second column, ch2
subplot(3,2,2); image(rewards_ch2_both, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch2); colormap('jet');    title('Ch 2');


subplot(3,2,4); image(rewards_ch2_470, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch2); colormap('jet');    


subplot(3,2,6); image(rewards_ch2_405, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_ch2); colormap('jet');    



if saveOn
    saveas(gcf, fullfile(savepath, 'random_rewards_alternateLEDs.fig'));
    saveas(gcf, fullfile(savepath, 'random_rewards_alternateLEDs.jpg'));
end


%% save plot of baseline fluorecense values across conditions
fbl = zeros(1,6);
fbl_sem = zeros(1,6);
fbl(1) = nanmean(TE.Photometry_ch1.data(1).blF_raw(LED12trials));
fbl(2) = nanmean(TE.Photometry_ch1.data(1).blF_raw(LED1trials));
fbl(3) = nanmean(TE.Photometry_ch1.data(1).blF_raw(LED2trials));
fbl(4) = nanmean(TE.Photometry_ch2.data(1).blF_raw(LED12trials));
fbl(5) = nanmean(TE.Photometry_ch2.data(1).blF_raw(LED1trials));
fbl(6) = nanmean(TE.Photometry_ch2.data(1).blF_raw(LED2trials));

fbl_sem(1) = nanSEM(TE.Photometry_ch1.data(1).blF_raw(LED12trials));
fbl_sem(2) = nanSEM(TE.Photometry_ch1.data(1).blF_raw(LED1trials));
fbl_sem(3) = nanSEM(TE.Photometry_ch1.data(1).blF_raw(LED2trials));
fbl_sem(4) = nanSEM(TE.Photometry_ch2.data(1).blF_raw(LED12trials));
fbl_sem(5) = nanSEM(TE.Photometry_ch2.data(1).blF_raw(LED1trials));
fbl_sem(6) = nanSEM(TE.Photometry_ch2.data(1).blF_raw(LED2trials));

% fbl(1) = nanmean(TE.Photometry.data(1).blF_raw(LED12trials));
% fbl(2) = nanmean(TE.Photometry.data(2).blF_raw(LED12trials));
% fbl(3) = nanmean(TE.Photometry.data(1).blF_raw(LED1trials));
% fbl(4) = nanmean(TE.Photometry_1alt.data(2).blF_raw(LED1trials));
% fbl(5) = nanmean(TE.Photometry_2alt.data(1).blF_raw(LED2trials));
% fbl(6) = nanmean(TE.Photometry.data(2).blF_raw(LED2trials));
% 
% fbl_sem(1) = nanSEM(TE.Photometry.data(1).blF_raw(LED12trials));
% fbl_sem(2) = nanSEM(TE.Photometry.data(2).blF_raw(LED12trials));
% fbl_sem(3) = nanSEM(TE.Photometry.data(1).blF_raw(LED1trials));
% fbl_sem(4) = nanSEM(TE.Photometry_1alt.data(2).blF_raw(LED1trials));
% fbl_sem(5) = nanSEM(TE.Photometry_2alt.data(1).blF_raw(LED2trials));
% fbl_sem(6) = nanSEM(TE.Photometry.data(2).blF_raw(LED2trials));
ensureFigure('baseline_fluor', 1);

bar(1:6, fbl, 'b'); hold on;
errorbar(1:6,fbl, fbl_sem, '.', 'LineWidth', 2 );
set(gca, 'XTick', 1:6, 'XTickLabel',...
    {'ch1 both', 'ch1 470', 'ch1 405', 'ch2 both', 'ch2 470', 'ch2 405'});
ylabel('Baseline Fluor (V)');
if saveOn
    saveas(gcf, fullfile(savepath, 'baseline_fluor.fig'));
    saveas(gcf, fullfile(savepath, 'baseline_fluor.jpg'));
end

%% given that you've used 