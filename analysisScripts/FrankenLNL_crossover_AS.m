saveOn = 1;
%% 
sessions = bpLoadSessions;

%%
TE = makeTE_FrankenLNL_4odors(sessions);

%%
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRCaMP(ch2)
channels=[]; dFFMode = {}; BL = {}; 
if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'expFit';
%     dFFMode{end+1} = 'simple';    
    BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
    BL{end + 1} = [1 4];    
end

%% process photometry
% baseline by trial
 TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue2', 'channels', channels, 'baseline', BL);
 
 %%
% basepath = uigetdir;
basepath = 'Z:\SummaryAnalyses\Franken_LNL_crossover';
sep = strfind(TE.filename{1}, '_');
if length(sessions) > 1
    subjectName = TE.filename{1}(1:sep(2)-1);
else
    subjectName = TE.filename{1}(1:end-4);
end
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

%%
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'trialType', 1));

%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%% cross sessions bleaching curve and exponential fits
for channel = channels
    figname = ['sessionBleach_Correction_ch' num2str(channel)];
    ensureFigure(figname, 1);
    plot(TE.Photometry.data(channel).blF_raw, 'k'); hold on;
    plot(TE.Photometry.data(channel).blF, 'r');
    if saveOn
        saveas(gcf, fullfile(savepath, [figname '.fig']));
        saveas(gcf, fullfile(savepath, [figname '.jpg']));
    end
    % cross trial bleaching fits for each session plotted as axis array
    if 1 %channel == 1
        figname = ['trialBleach_Correction_ch' num2str(channel)];
        f1 = ensureFigure(figname, 1);
        nSessions = size(TE.Photometry.bleachFit, 1);
        subA = ceil(sqrt(nSessions));
        for counter = 1:nSessions
            subplot(subA, subA, counter);
            plot(TE.Photometry.bleachFit(counter, channel).trialTemplateFullX, TE.Photometry.bleachFit(counter, channel).trialTemplateFull, 'g'); hold on;            
            plot(TE.Photometry.bleachFit(counter, channel).trialTemplateX, TE.Photometry.bleachFit(counter, channel).trialTemplate, 'b'); 
            plot(TE.Photometry.bleachFit(counter, channel).fitX, TE.Photometry.bleachFit(counter, channel).trialFit, 'r');
        %     title(num2str(counter));    
        end
        % average of all trials for this channel to eyeball correction
        figname2 = ['corrected_allTrials_ch' num2str(channel)]; 
        ensureFigure(figname2, 1);
        [ha, hl] = phPlotAverageFromTE(TE, 1:length(TE.filename), channel,...
    'FluorDataField', 'ZS', 'window', [0.1, max(TE.Photometry.xData) - min(TE.Photometry.xData)], 'zeroTimes', TE.Photometry.startTime); %high value, reward
        if saveOn
            saveas(f1, fullfile(savepath, figname), 'fig');
            saveas(f1, fullfile(savepath, figname), 'jpeg');            
            saveas(gcf, fullfile(savepath, [figname2 '.fig']));
            saveas(gcf, fullfile(savepath, [figname2 '.jpg']));            
        end
    end
end

%% 
frankenLNL_conditions;

%% raster plots
PhotometryField = 'Photometry';
CLimFactor = 3;
for channel = channels
    saveName = ['Crossover_phRasters_ch' num2str(channel)];
    ensureFigure(saveName, 1);
    
    subplot(1,3,1);
    phRasterFromTE(TE, trialsByType{1}, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-4 7], 'zeroTimes', TE.Cue1);
    xlabel('time from cue (s)');  title('A->B');
    subplot(1,3,2);
    phRasterFromTE(TE, trialsByType{2}, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-4 7], 'zeroTimes', TE.Cue1);
    xlabel('time from cue (s)');  title('C->D');    
    subplot(1,3,3);
    phRasterFromTE(TE, trialsByType{3} & rewardTrials, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-4 7], 'zeroTimes', TE.Cue1);
    xlabel('time from cue (s)');  title('uncued reward?');       
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end
end

%% lick raster



    saveName = 'Crossover_lickRaster';
    ensureFigure(saveName, 1);
    
    subplot(1,3,1);
%     
%         eventRasterFromTE(TE, csPlusTrials, 'Port1In', 'trialNumbering', 'global',...
%         'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    
    window = [-4 7];
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', 'consecutive', 'window', window, 'zeroTimes', TE.Cue1);
    xlabel('time from cue (s)');  title('A->B'); set(gca, 'XLim', window);
    subplot(1,3,2);
    eventRasterFromTE(TE, trialsByType{2}, 'Port1In',  'trialNumbering', 'consecutive', 'window', window, 'zeroTimes', TE.Cue1);
    xlabel('time from cue (s)');  title('C->D');    set(gca, 'XLim', window);
    subplot(1,3,3);
    eventRasterFromTE(TE, trialsByType{3} & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive', 'window', window, 'zeroTimes', TE.Cue1);
    xlabel('time from cue (s)');  title('uncued reward?');       set(gca, 'XLim', window);
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end
