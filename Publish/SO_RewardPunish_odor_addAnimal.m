% SO_RewardPunish_odor_addAnnimal

DB = dbLoadExperiment('SO_RewardPunish_odor');
saveOn = 1;
%%
sessions = bpLoadSessions;
%%
TE = makeTE_SO_RewardPunish_Odor_v2(sessions);

%%
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRCaMP(ch2)
channels=[]; dFFMode = {}; BL = {};

if sessions(1).SessionData.TrialSettings(1).GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'expFit';
%     dFFMode{end+1} = 'simple';    
    BL{end + 1} = [0 sessions(1).SessionData.TrialSettings(1).PreCsRecording];
end

if sessions(1).SessionData.TrialSettings(1).GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
    BL{end + 1} = [0 sessions(1).SessionData.TrialSettings(1).PreCsRecording];
end

%% process photometry
% baseline by trial
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL);
% baseline expfit
TE.PhotometryExpFit = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'expFit', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL);

%%
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
DB = dbRegisterAnimal(DB, subjectName);
savepath = dbGetAnimalPath(DB, subjectName);

%% truncate sessions
% usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.Omi, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.Delay, 'referenceFromEnd', 1);
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'ReinforcementOutcome', 'Reward'));
%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% cross sessions bleaching curve and exponential fits
for channel = channels
    figname = ['sessionBleach_Correction_ch' num2str(channel)];
    ensureFigure(figname, 1);
    plot(TE.PhotometryExpFit.data(channel).blF_raw, 'k'); hold on;
    plot(TE.PhotometryExpFit.data(channel).blF, 'r');
    if saveOn
        saveas(gcf, fullfile(savepath, [figname '.fig']));
        saveas(gcf, fullfile(savepath, [figname '.jpg']));
    end
    % cross trial bleaching fits for each session plotted as axis array
        figname = ['trialBleach_Correction_ch' num2str(channel)];
        ensureFigure(figname, 1);
        nSessions = size(TE.PhotometryExpFit.bleachFit, 1);
        subA = ceil(sqrt(nSessions));
        for counter = 1:nSessions
            subplot(subA, subA, counter);
            plot(TE.PhotometryExpFit.bleachFit(counter, channel).trialTemplateFullX, TE.PhotometryExpFit.bleachFit(counter, channel).trialTemplateFull, 'g'); hold on;            
            plot(TE.PhotometryExpFit.bleachFit(counter, channel).trialTemplateX, TE.PhotometryExpFit.bleachFit(counter, channel).trialTemplate, 'b'); 
            plot(TE.PhotometryExpFit.bleachFit(counter, channel).fitX, TE.PhotometryExpFit.bleachFit(counter, channel).trialFit, 'r');
        %     title(num2str(counter));    
        end
        % average of all trials for this channel to eyeball correction
        figname2 = ['corrected_allTrials_ch' num2str(channel)]; 
        ensureFigure(figname2, 1);
        [ha, hl] = phPlotAverageFromTE(TE, 1:length(TE.filename), channel,...
    'FluorDataField', 'ZS', 'window', [0.1, max(TE.PhotometryExpFit.xData) - min(TE.PhotometryExpFit.xData)], 'zeroTimes', TE.PhotometryExpFit.startTime); %high value, reward
        if saveOn
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.fig'));
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.jpg'));
            saveas(gcf, fullfile(savepath, [figname2 '.fig']));
            saveas(gcf, fullfile(savepath, [figname2 '.jpg']));            
        end
end
