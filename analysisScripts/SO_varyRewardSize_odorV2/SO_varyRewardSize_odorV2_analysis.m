saveOn = 1;
sessions = bpLoadSessions;
%%
TE = makeTE_SO_varyRewardSize_odorV2(sessions);

%%
photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
climfactor = 2;
window = [-6 4];

ensureDirectory(savepath);
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
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'DeliverStimulus', 'channels', channels, 'baseline', BL);
% baseline expfit
% TE.PhotometryExpFit = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'expFit', 'zeroField', 'DeliverStimulus', 'channels', channels, 'baseline', BL);

%% truncate sessions
% usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.Omi, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.Delay, 'referenceFromEnd', 1);
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'ReinforcementOutcome', 'Reward'));
%% animal name and path


basepath = uigetdir;
% basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\';
% basepath = 'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\';
% basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_20and17_combined\';
sep = strfind(TE.filename{1}, '_');
animal = TE.filename{1}(1:sep(2)-1);
disp(animal);
savepath = fullfile(basepath, animal);
ensureDirectory(savepath);
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% trial types
SO_varyRewardSize_odorV2_conditions;

%% photometry averages, zscored
%     ylim = [-2 8];
    window = [-4 2];
    channel = 1;
    fdField = 'ZS';
    saveName = sprintf('Averages_%s', animal);  
    ensureFigure(saveName, 1);
    
    subplot(2,1,1); hold on;
    varargin = {'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'b', 'c', 'k'}};
    axh = [];
    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1:3]), 'Port1In', varargin{:});
    varargin = {'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'b--', 'c--', 'k--'}};
    axh = [];
    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([4:6]), 'Port1In', varargin{:});    
%     legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    ylabel('licks (1/s)'); xlabel('time from reinforcement (s)'); 
    
    subplot(2, 1, 2); hold on;
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType(1:3), channel,...
        'PhotometryField', 'Photometry', 'FluorDataField', fdField, 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'b', 'c','k'});
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType(4:6), channel,...
        'PhotometryField', 'Photometry', 'FluorDataField', fdField, 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'b--', 'c--','k--'});
%     legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    title('Reward'); ylabel(sprintf('BF %s', fdField)); textBox(animal, gca, [0.25 0.95], 8);%set(gca, 'YLim', ylim);
    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end
    
    %% raster plots
    saveName = sprintf('rasters_%s', animal); 
    ensureFigure(saveName, 1);
    mcLandscapeFigSetup(gcf); 

    trialTypes = [1 2 5 3 4 6];
    window = [-6 3];
    titles  = {'cued small', 'cued big', 'omit', 'uncued small', 'uncued big', 'null'};
    for typeCounter = 1:length(trialTypes)
        trialType = trialTypes(typeCounter);
        sessionChanges = find(diff(TE.sessionIndex(trialsByType{trialType}, :))) + 1;
        subplot(2,6,typeCounter);
        phRasterFromTE(TE, trialsByType{trialType}, channel, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor,...
            'zeroTimes', TE.usZeros, 'window', window, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,       
        title(titles{typeCounter}); 
%         cellfun(@(x) x(1), TE.Cue)
        
        subplot(2,6,typeCounter + 6);
        eventRasterFromTE(TE, trialsByType{trialType}, 'Port1In', 'trialNumbering', 'consecutive',...
            'zeroTimes', TE.usZeros, 'window', window);
        line(repmat([window(1); window(2)], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
        set(gca, 'XLim', window); 
    end
    subplot(2,6,8); xlabel('Time from Reinforcement (s)');
    subplot(2,6,1); t = textBox(animal, gca, [0.1 0.95]); set(t, 'Color', [1 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end