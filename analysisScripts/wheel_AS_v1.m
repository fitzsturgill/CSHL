% 4/10/17  Analysis script for pavlovian reversals using LickNoLick_Odor_V2
% protocol



saveOn = 0;
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
    dFFMode{end+1} = 'expFit';
%     BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
%     BL{end + 1} = [2 4];    
end


TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'expFit',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [0 2], 'startField', 'Baseline');

%%
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', 30, 'Fs', 20, 'startField', 'Start');
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
%%