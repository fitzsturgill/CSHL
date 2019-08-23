% wheel_addAnimal

DB = dbLoadExperiment('wheel');
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

try
    baselineEnd = sessions(1).SessionData.NidaqData{1,2}.duration(1) - 0.2;
catch
    baselineEnd = (size(sessions.SessionData.NidaqData{1,1}, 1) / 6100) - 0.2;
end
        
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [0 baselineEnd], 'startField', 'Baseline', 'downsample', 305);

%%
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', baselineEnd + 1, 'Fs', 20, 'startField', 'Start');
%% pupil data
%  [wheelY_new, wheelTimes_new] = resample(wheelY, wheelTimes, 20, 'linear');
folderSuffix = ''; % or enter folder suffix on command line
%%
TE = addPupilometryToTE(TE, 'duration', baselineEnd + 1, 'zeroField', 'Baseline', 'startField', 'Baseline', 'frameRate', 60, 'frameRateNew', 20, 'folderSuffix', folderSuffix, 'fillNaNs', 1);
%%
TE.Whisk = addWhiskingToTE(TE, 'duration', baselineEnd + 1, 'zeroField', 'Baseline', 'startField', 'Baseline', 'sampleRate', 60, 'sampleRateNew', 20, 'folderSuffix', folderSuffix);

%%
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
DB = dbRegisterAnimal(DB, subjectName);
savepath = dbGetAnimalPath(DB, subjectName);

%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%%
figure; subplot(2,2,1); imagesc(TE.Photometry.data(1).raw); subplot(2,2,2); imagesc(TE.Photometry.data(1).ZS); subplot(2,2,3); imagesc(TE.Wheel.data.V); subplot(2,2,4); imagesc(TE.pupil.pupDiameterNorm);
