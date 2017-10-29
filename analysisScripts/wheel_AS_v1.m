% wheel analysis script

saveOn = 0;
%%
saveOn = 1;
%%
sessions = bpLoadSessions; % load sessions
%% 
TE = makeTE_wheel_v1(sessions); % make TE structure
%%
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRGECO1a (ch2)
channels=[]; dFFMode = {}; %BL = {};
if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'simple';
%     BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
%     BL{end + 1} = [2 4];    
end

baselineEnd = 119;
% baselineEnd = 29;
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [0 baselineEnd], 'startField', 'Baseline', 'downsample', 305);

%%
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', baselineEnd + 1, 'Fs', 20, 'startField', 'Start');

%% pupil data
%  [wheelY_new, wheelTimes_new] = resample(wheelY, wheelTimes, 20, 'linear');
folderSuffix = ''; % or enter folder suffix on command line
%%
TE = addPupilometryToTE(TE, 'duration', baselineEnd + 1, 'zeroField', 'Baseline', 'startField', 'Baseline', 'frameRate', 60, 'frameRateNew', 20, 'folderSuffix', folderSuffix);
%% Now saved in directory according to first session filename
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
% basepath = uigetdir;

basepath = 'Z:\SummaryAnalyses\wheel_v1\';
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

%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%% plot raw and smoothed scatter plots of all the data (excepting the first few trials)
nPoints = numel(TE.Photometry.data(2).ZS(1:end,:)); 
ensureFigure('scatter', 1); 
ChAT_raw = TE.Photometry.data(1).ZS(:);
DAT_raw = TE.Photometry.data(2).ZS(:);
a=zeros(2,1);
smoothfactor = 100;
a(1) = subplot(1,2,1); scatter(DAT_raw, ChAT_raw, '.'); xlabel('DAT fluor (Zscored)'); ylabel('ChAT fluor (Zscored)'); title('raw');
a(2) = subplot(1,2,2); scatter(smooth(DAT_raw, smoothfactor), smooth(ChAT_raw, smoothfactor), '.'); xlabel('DAT fluor (Zscored)'); ylabel('ChAT fluor (Zscored)'); title('smoothed');
sameXYScale(a); %sameXScale(a);sameYScale(a);
% setXYsameLimit(a(1), 0);setXYsameLimit(a(2), 0);

if saveOn
    saveas(gcf, fullfile(savepath, 'scatter.fig'));
    saveas(gcf, fullfile(savepath, 'scatter.jpg'));
end
%%
window = [-2 2];
fs = 20; % sample rate
blSamples = (0 - window(1)) * fs;
[rewards_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2]);
rewards_chat = extractDataByTimeStamps(TE.Photometry.data(1).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2]);

% local dFF
bl_dat = nanmean(rewards_dat(:,1:blSamples), 2);
rewards_dat = bsxfun(@minus, rewards_dat, bl_dat);
rewards_dat = bsxfun(@rdivide, rewards_dat, bl_dat);
sd_dat = nanmean(nanstd(rewards_dat(:,1:blSamples)));
bl_chat = nanmean(rewards_chat(:,1:blSamples), 2);
rewards_chat = bsxfun(@minus, rewards_chat, bl_chat);
rewards_chat = bsxfun(@rdivide, rewards_chat, bl_chat);
sd_chat = nanmean(nanstd(rewards_chat(:,1:blSamples)));
nTrials = length(TE.filename);
ts_abs = zeros(size(ts));
for counter = 1:length(ts)
    ts_abs(counter) = ts(counter) + TE.trialStartTimeStamp(tn(counter));    
end

iri_pre = [Inf; diff(ts_abs)];
iri_post = [diff(ts_abs); Inf];

[~, I] = sort(iri_pre);
iri_pre_sorted = iri_pre(I);
rewards_dat_sorted = rewards_dat(I, :);
rewards_chat_sorted = rewards_chat(I, :);
climFactor = 4;
clim_chat = [-climFactor * sd_chat, climFactor * sd_chat];% + 0.003;
clim_dat = [-climFactor * sd_dat, climFactor * sd_dat];% + 0.03;
ensureFigure('random_rewards', 1); 
    % subplot(3,2,1); image(rewards_chat, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); colormap('jet');  title('ChAT');
% subplot(3,2,2); image(rewards_dat, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', [min(min(rewards_dat_sorted)), max(max(rewards_dat_sorted))]); colormap('jet');    title('DAT');
subplot(3,2,1); image(rewards_chat, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_chat); colormap('jet');  title('ChAT');
subplot(3,2,2); image(rewards_dat, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_dat); colormap('jet');    title('DAT');
xdata = linspace(window(1), window(2), size(rewards_chat, 2));
subplot(3,2,3); plot(xdata, nanmean(rewards_chat));
subplot(3,2,4); plot(xdata, nanmean(rewards_dat));
subplot(3,2,5); triggeredEventRasterFromTE(TE, 'Port1In', TE.Reward);
set(gca, 'YLim', [0 size(rewards_chat, 1)]);

% % reward licks vs. trial number to truncate
% rl = extractTriggeredEvents(TE, 'Port1In', TE.Reward);
% ensureFigure('truncate', 1);
% rl_trials = unique(rl.eventTrials);
% rl_count = zeros(size(rl_trials));
% for counter = 1:length(rl_trials)
%     trial = rl_trials(counter);
%     rl_count(counter) = sum(rl.eventTimes > 0 & rl.eventTrials == trial);    
% end
% plot(rl_trials, smooth(rl_count)); ylabel('# reward licks'); xlabel('trial #');
    

if saveOn
    saveas(gcf, fullfile(savepath, 'random_rewards.fig'));
    saveas(gcf, fullfile(savepath, 'random_rewards.jpg'));
end


%%
% good trials for 5/31 session,  1, 2, 3, 17
% good trials for 6/1 session 1, 4-  !!!! trial 1 is good for showing DAT
% and ChAT correlations with reward and without but doesn't have nice pupil
% diameter
% good trials with pupil traces that needed gap filling: 12
trial = 5; % 7;
ensureFigure('examples', 1);
subplot(4,1,1);
ydata = TE.Photometry.data(1).raw(trial, :);    
plot(TE.Photometry.xData, ydata, 'k'); hold on;
tsx = repmat(TE.Reward{trial}(:,1), 1, 2)';
tsy = [repmat(min(min(ydata)), 1, size(tsx, 2)); repmat(max(max(ydata)), 1, size(tsx, 2))];    
plot(tsx, tsy, 'r'); ylabel('ChAT');
subplot(4,1,2);
ydata = TE.Photometry.data(2).raw(trial, :);    
plot(TE.Photometry.xData, ydata, 'k'); hold on;
tsx = repmat(TE.Reward{trial}(:,1), 1, 2)';
tsy = [repmat(min(min(ydata)), 1, size(tsx, 2)); repmat(max(max(ydata)), 1, size(tsx, 2))];    
plot(tsx, tsy, 'r'); ylabel('DAT');
pupField = 'pupDiameter';
try
    subplot(4,1,3); plot(TE.pupil.xData, TE.pupil.(pupField)(trial, :)); ylabel('Pupil Diameter');
    set(gca, 'YLim', [percentile(TE.pupil.(pupField)(trial, :), 0.03), percentile(TE.pupil.(pupField)(trial, :), 0.97)]);
catch
end
subplot(4,1,4); plot(TE.Wheel.xData, TE.Wheel.data.V(trial, :)); ylabel('Velocity');

if saveOn
    saveas(gcf, fullfile(savepath, 'examples.fig'));
    saveas(gcf, fullfile(savepath, 'examples.jpg'));
end

%% McKnight special for 6/1 session, 
% good trials for 5/31 session,  1, 2, 3, 17
% good trials for 6/1 session 1, 4-  !!!! trial 1 is good for showing DAT
% and ChAT correlations with reward and without but doesn't have nice pupil
% diameter
% good trials with pupil traces that needed gap filling: 12
trial = 6;
ensureFigure('examples', 1);
subplot(2,1,1);
ydata = TE.Photometry.data(1).ZS(trial, :);    
plot(TE.Photometry.xData, ydata, 'g'); hold on;
tsx = repmat(TE.Reward{trial}(:,1), 1, 2)';
tsy = [repmat(min(min(ydata)), 1, size(tsx, 2)); repmat(max(max(ydata)), 1, size(tsx, 2))];    
plot(tsx, tsy, 'b'); ylabel('ChAT (ZScore)');
subplot(2,1,2);
ydata = TE.Photometry.data(2).ZS(trial, :);    
plot(TE.Photometry.xData, ydata, 'r'); hold on;
tsx = repmat(TE.Reward{trial}(:,1), 1, 2)';
tsy = [repmat(min(min(ydata)), 1, size(tsx, 2)); repmat(max(max(ydata)), 1, size(tsx, 2))];    
plot(tsx, tsy, 'b'); ylabel('DAT (ZScore)');
if saveOn
    saveas(gcf, fullfile(savepath, 'ChAT_vs_DAT_randomReward_example.fig'));
    saveas(gcf, fullfile(savepath, 'ChAT_vs_DAT_randomReward_example.jpg'));
    saveas(gcf, fullfile(savepath, 'ChAT_vs_DAT_randomReward_example.epsc'));
end

%% phase analysis (hilbert transform)
trial = 1;
trials = 1:length(TE.filename);
Fs = 20;
bp = [0.1 2];


TE.timeFromReward = bpCalcTimeFromEvent(TE, 'Reward', 'dataStart', TE.Photometry.startTime, 'trialStart', TE.trialStartTimeStamp, 'duration', baselineEnd + 1);

% de-trend the signal with a band-pass filter
% you may need the signal processing toolbox ....
bp = bp * 2 / Fs; % convert Hz to radians/S
[N, Wn] = buttord( bp, bp .* [.5 1.5], 3, 20); 
[B,A] = butter(N,Wn);
sizeData = size(TE.Photometry.data(1).raw);
nSamples = numel(TE.Photometry.data(1).raw);
hdata = struct(...
    'data', zeros(sizeData),...
    'filtData', zeros(sizeData),...
    'hilb', zeros(sizeData),...
    'phase', zeros(sizeData),...
    'amp', zeros(sizeData)...
    );
hdata = repmat(hdata, size(channels));
for channel = 1:2
    for trial = 1:length(TE.filename)
        hdata(channel).data(trial,:) = TE.Photometry.data(channel).ZS(trial, :);    
        hdata(channel).filtData(trial, :) = filtfilt(B,A,hdata(channel).data(trial,:)); % zero-phase filtering
        hdata(channel).hilb(trial, :) = hilbert(hdata(channel).filtData(trial,:));
        hdata(channel).phase(trial, :) = angle(hdata(channel).hilb(trial,:));
        hdata(channel).amp(trial, :) = abs(hdata(channel).hilb(trial,:));
    end
end
pm = [4 1];
trial = 1;
ensureFigure('Hilbert_examples', 1);
subplot(pm(1), pm(2), 1); plot(hdata(1).data(trial,:), 'g'); hold on; plot(hdata(2).data(trial,:), 'r');
subplot(pm(1), pm(2), 2); plot(hdata(1).filtData(trial,:), 'g'); plot(hdata(2).filtData(trial,:), 'r');
subplot(pm(1), pm(2), 3); plot(hdata(1).phase(trial,:), 'g'); hold on; plot(hdata(2).phase(trial,:), 'r');
subplot(pm(1), pm(2), 4); plot(hdata(1).amp(trial,:), 'g'); hold on; plot(hdata(2).amp(trial,:), 'r');

ensureFigure('Hilbert_scatter', 1);
% subplot(2,2,1); scatter(reshape(TE.timeFromReward, nSamples, 1), reshape(hdata(1).phase, nSamples, 1), '.', 'MarkerFaceColor', 'g'); 
% subplot(2,2,2); scatter(reshape(TE.timeFromReward, nSamples, 1), reshape(hdata(2).phase, nSamples, 1), '.', 'MarkerFaceColor', 'r');
% subplot(2,2,3); scatter(reshape(TE.timeFromReward, nSamples, 1), reshape(hdata(1).amp, nSamples, 1), '.', 'MarkerFaceColor', 'g'); 
% subplot(2,2,4); scatter(reshape(TE.timeFromReward, nSamples, 1), reshape(hdata(2).amp, nSamples, 1), '.', 'MarkerFaceColor', 'r');
% subplot(1,2,1); scatter3(reshape(TE.timeFromReward, nSamples, 1), reshape(hdata(1).phase, nSamples, 1), reshape(hdata(2).phase, nSamples, 1)); %, '.', 'MarkerFaceColor', 'r');
% subplot(1,2,2); scatter3(reshape(TE.timeFromReward, nSamples, 1), reshape(hdata(1).amp, nSamples, 1), reshape(hdata(2).amp, nSamples, 1)); %, '.', 'MarkerFaceColor', 'g'); 
subplot(1,2,1);
scatter(TE.timeFromReward(:), hdata(1).phase(:) - hdata(2).phase(:), '.');
subplot(1,2,2);
histogram2(TE.timeFromReward, hdata(1).phase - hdata(2).phase);

ampPercentileForThresh = 0.2;
ampThresh1 = percentile(hdata(1).amp(trials, :), ampPercentileForThresh);
ampThresh2 = percentile(hdata(2).amp(trials, :), ampPercentileForThresh);

timeFromReward = TE.timeFromReward(trials, :);
phaseDiff = hdata(1).phase(trials, :) - hdata(2).phase(trials, :);
ensureFigure('binnedPhases', 1);
validPoints = ~isinf(timeFromReward) & hdata(1).amp(trials, :) > ampThresh1 & hdata(2).amp(trials, :) > ampThresh2;
[phaseMeans, phaseErrors, timeFromReward] = binnedMeansXY(timeFromReward(validPoints), phaseDiff(validPoints), 20);
errorbar(timeFromReward, phaseMeans, phaseErrors); xlabel('time from reward'); ylabel('hilbert phase');
if saveOn
    saveas(gcf, fullfile(savepath, 'binnedPhases.fig'));
    saveas(gcf, fullfile(savepath, 'binnedPhases.jpg'));
    saveas(gcf, fullfile(savepath, 'binnedPhases.epsc'));
end