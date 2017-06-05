% 4/10/17  Analysis script for pavlovian reversals using LickNoLick_Odor_V2
% protocol

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


TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'bySession',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [0 2], 'startField', 'Baseline', 'downsample', 305);

%%
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', 30, 'Fs', 20, 'startField', 'Start');

%% pupil data
%  [wheelY_new, wheelTimes_new] = resample(wheelY, wheelTimes, 20, 'linear');

TE = addPupilometryToTE(TE, 'duration', 30, 'zeroField', 'Baseline', 'startField', 'Baseline', 'frameRate', 60, 'frameRateNew', 20);
%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = uigetdir;
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

%% plot raw and smoothed scatter plots of all the data (excepting the first few trials)
nPoints = numel(TE.Photometry.data(2).ZS(4:end,:)); 
ensureFigure('scatter', 1); 
ChAT_raw = reshape(TE.Photometry.data(1).ZS(4:end,:), nPoints, 1);
DAT_raw = reshape(TE.Photometry.data(2).ZS(4:end,:), nPoints, 1);
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
% good trials,  7
% good trials with pupil traces that needed gap filling: 12
trial = 7;
ensureFigure('examles', 1);
subplot(4,1,1); plot(TE.Photometry.xData, TE.Photometry.data(1).raw(trial,:)); ylabel('ChAT');
subplot(4,1,2); plot(TE.Photometry.xData, TE.Photometry.data(2).raw(trial,:)); ylabel('DAT');
subplot(4,1,3); plot(TE.pupil.xData, TE.pupil.pupDiameter(trial, :)); ylabel('Pupil Diameter');
set(gca, 'YLim', [percentile(TE.pupil.pupDiameter(trial, :), 0.03), percentile(TE.pupil.pupDiameter(trial, :), 0.97)]);
subplot(4,1,4); plot(TE.Wheel.xData, TE.Wheel.data.V(trial, :)); ylabel('Velocity');

if saveOn
    saveas(gcf, fullfile(savepath, 'examples.fig'));
    saveas(gcf, fullfile(savepath, 'examples.jpg'));
end