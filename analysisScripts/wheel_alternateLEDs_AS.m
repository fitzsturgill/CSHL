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
basepath = uigetdir;
sep = strfind(TE.filename{1}, '.');
subjectName = TE.filename{1}(1:sep(1)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

%%
channels = [1 2];
dFFMode = {'simple', 'simple'};

TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [0 29], 'startField', 'Baseline', 'downsample', 305);

% channel 1 demodulated via channel 2 reference
TE.Photometry_1alt = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [0 29], 'startField', 'Baseline', 'downsample', 305, 'refChannels', [1 1]);

% channel 2 demodulated via channel 1 reference
TE.Photometry_2alt = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [0 29], 'startField', 'Baseline', 'downsample', 305, 'refChannels', [2 2]);

%%

LED1trials = TE.LED1_amp > 0 & TE.LED2_amp == 0;
LED2trials = TE.LED2_amp > 0 & TE.LED1_amp == 0;
LED12trials = TE.LED1_amp > 0 & TE.LED2_amp > 0;

%%

window = [-2 2];
fs = 20; % sample rate
blSamples = (0 - window(1)) * fs;

% LED 1 and 2
[rewards_dat_12, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2], LED12trials);
rewards_chat_12 = extractDataByTimeStamps(TE.Photometry.data(1).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2], LED12trials);
bl_dat_12 = nanmean(rewards_dat_12(:,1:blSamples), 2);
rewards_dat_12 = bsxfun(@minus, rewards_dat_12, bl_dat_12);
bl_chat_12 = nanmean(rewards_chat_12(:,1:blSamples), 2);
rewards_chat_12 = bsxfun(@minus, rewards_chat_12, bl_chat_12);

% LED 1 only
[rewards_dat_1, ts, tn] = extractDataByTimeStamps(TE.Photometry_1alt.data(2).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2], LED1trials);
rewards_chat_1 = extractDataByTimeStamps(TE.Photometry_1alt.data(1).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2], LED1trials);
bl_dat_1 = nanmean(rewards_dat_1(:,1:blSamples), 2);
rewards_dat_1 = bsxfun(@minus, rewards_dat_1, bl_dat_1);
bl_chat_1 = nanmean(rewards_chat_1(:,1:blSamples), 2);
rewards_chat_1 = bsxfun(@minus, rewards_chat_1, bl_chat_1);

% LED 2 only
[rewards_dat_2, ts, tn] = extractDataByTimeStamps(TE.Photometry_2alt.data(2).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2], LED2trials);
rewards_chat_2 = extractDataByTimeStamps(TE.Photometry_2alt.data(1).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2], LED2trials);
bl_dat_2 = nanmean(rewards_dat_2(:,1:blSamples), 2);
rewards_dat_2 = bsxfun(@minus, rewards_dat_2, bl_dat_2);
bl_chat_2 = nanmean(rewards_chat_2(:,1:blSamples), 2);
rewards_chat_2 = bsxfun(@minus, rewards_chat_2, bl_chat_2);

% calc mean and SD for LED 1 and 2 trials
sd_chat = nanmean(nanstd(rewards_chat_12(:,1:blSamples)));
sd_dat = nanmean(nanstd(rewards_dat_12(:,1:blSamples)));

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
clim_chat = [-climFactor * sd_chat, climFactor * sd_chat];
clim_dat = [-climFactor * sd_dat, climFactor * sd_dat];

ensureFigure('random_rewards_alternateLEDs', 1); 
subplot(3,2,1); image(rewards_chat_12, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_chat); colormap('jet');  title('ChAT');
subplot(3,2,3); image(rewards_chat_1, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_chat); colormap('jet');  ylabel('LED1 only');
subplot(3,2,5); image(rewards_chat_2, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_chat); colormap('jet');  ylabel('LED2 only');

subplot(3,2,2); image(rewards_dat_12, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_dat); colormap('jet');    title('DAT');
subplot(3,2,4); image(rewards_dat_1, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_dat); colormap('jet');    
subplot(3,2,6); image(rewards_dat_2, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', clim_dat); colormap('jet');    


if saveOn
    saveas(gcf, fullfile(savepath, 'random_rewards_alternateLEDs.fig'));
    saveas(gcf, fullfile(savepath, 'random_rewards_alternateLEDs.jpg'));
end
