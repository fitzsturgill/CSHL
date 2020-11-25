% DAT_acquisition_special
saveOn = 1;
filenames =     {...    
    'DAT_2_SO_varyRewardSize_odorV2_Apr01_2016_Session1.mat',...
    'DAT_2_SO_varyRewardSize_odorV2_Apr02_2016_Session1.mat',...
    'DAT_2_SO_RewardPunish_odor_Apr03_2016_Session1.mat',...
    'DAT_2_SO_RewardPunish_odor_Apr03_2016_Session2.mat',...
    'DAT_2_SO_RewardPunish_odor_Apr04_2016_Session2.mat',...
    'DAT_2_SO_RewardPunish_odor_Apr05_2016_Session1.mat'...
    };
filepaths = vertcat(...
 repmat({'Z:\FitzRig1\Data\DAT_2\SO_varyRewardSize_odorV2\Session Data\'}, 2, 1),...
 repmat({'Z:\FitzRig1\Data\DAT_2\SO_RewardPunish_odor\Session Data\'}, 4, 1)...
 );
sessions = bpLoadSessions([], filenames, filepaths);

%%
TE = makeTE_DAT2_special(sessions);
%%
channels=[]; dFFMode = {}; BL = {};

if sessions(1).SessionData.TrialSettings(1).GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'simple';
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


%%
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.Delay, 'referenceFromEnd', 1);
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'ReinforcementOutcome', 'Reward'));

%%
savepath = 'Z:\SummaryAnalyses\ChAT_PE_Manuscript\SO_RewardPunish_odor\DAT2_special\';
animal = 'DAT_2';
disp(animal);
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% load TE if it's already made
savepath = 'Z:\SummaryAnalyses\ChAT_PE_Manuscript\SO_RewardPunish_odor\DAT2_special\';
animal = 'DAT_2';
disp(animal);

load(fullfile(savepath, 'TE.mat'));
disp(['*** Loaded: ' fullfile(savepath, 'TE.mat')]);

%% kludge conditions
csPlusTrials = TE.odorValve == 5 & TE.reject == 0;
csMinusTrials = TE.odorValve == 6  & TE.reject == 0;
uncuedTrials = TE.odorValve == 0 & TE.reject == 0;
smallRewardTrials = TE.rewardSize == 2 & ismember(TE.ReinforcementOutcome, 'Reward')' & TE.reject == 0;
bigRewardTrials = TE.rewardSize == 8 & ismember(TE.ReinforcementOutcome, 'Reward')' & TE.reject == 0;
omitTrials = ismember(TE.ReinforcementOutcome, 'Omit')' & TE.reject == 0;

%% photometry averages, zscored
%     ylim = [-2 8];
    window = [-4 2];
    channel = 1;
    fdField = 'ZS';
    saveName = sprintf('Averages_%s', animal);  
    ensureFigure(saveName, 1);
    trialsToPlot = {...
        csPlusTrials & bigRewardTrials & ismember(TE.sessionIndex, 3:6);...
        uncuedTrials & bigRewardTrials & ismember(TE.sessionIndex, 3:6);...
        };
    subplot(2,1,1); hold on;
    varargin = {'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'b', 'c'}};
    axh = [];
    [ha, hl] = plotEventAverageFromTE(TE, trialsToPlot, 'Port1In', varargin{:});
    
%     legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    ylabel('licks (1/s)'); xlabel('time from reinforcement (s)'); 
    
    subplot(2, 1, 2); hold on;
    [ha, hl] = phPlotAverageFromTE(TE, trialsToPlot, channel,...
        'PhotometryField', 'Photometry', 'FluorDataField', fdField, 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'b', 'c'});
%     legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    title('Reward'); ylabel(sprintf('BF %s', fdField)); textBox(animal, gca, [0.25 0.95], 8);%set(gca, 'YLim', ylim);
    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end

%%
climfactor = 3;  
fdField = 'ZS';
linewidth = 4;
window = [-3 3];
offset = -0.1;
photometryField = 'Photometry';
    saveName = sprintf('rasters_%s', animal); 
    ensureFigure(saveName, 1);
    mcLandscapeFigSetup(gcf); 
%     trialsToPlot = repmat({1:length(TE.filename)}, 2, 1);
    trialsToPlot = {...
        csPlusTrials & bigRewardTrials;...
        uncuedTrials & bigRewardTrials;...
        };    
    window = [-6 3];
    titles  = {'cued big', 'uncued big', 'omit'};
    for typeCounter = 1:length(trialsToPlot)
        trials = trialsToPlot{typeCounter};
        sessionChanges = find(diff(TE.sessionIndex(trials, :))) + 1;
        subplot(2,3,typeCounter);
        phRasterFromTE(TE, trials, channel, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor,...
            'zeroTimes', TE.usZeros, 'window', window, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,       
        title(titles{typeCounter}); 
%         cellfun(@(x) x(1), TE.Cue)
        
        subplot(2,3,typeCounter + 3);
        eventRasterFromTE(TE, trials, 'Port1In', 'trialNumbering', 'consecutive',...
            'zeroTimes', TE.usZeros, 'window', window);
        line(repmat([window(1); window(2)], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
        set(gca, 'XLim', window); 
    end
    subplot(2,3,5); xlabel('Time from Reinforcement (s)');
    subplot(2,3,4); t = textBox(animal, gca, [0.1 0.95]); set(t, 'Color', [1 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end


%% for figure, like GCaMP
earlySessions = {'DAT_2_SO_varyRewardSize_odorV2_Apr01_2016_Session1.mat'};
midSessions = {'DAT_2_SO_varyRewardSize_odorV2_Apr02_2016_Session1.mat'};
lateSessions = {'DAT_2_SO_RewardPunish_odor_Apr04_2016_Session2.mat', 'DAT_2_SO_RewardPunish_odor_Apr05_2016_Session1.mat'};

earlyTrials = csPlusTrials & bigRewardTrials & TE.reject == 0 & ismember(TE.filename, earlySessions);
midTrials = csPlusTrials & bigRewardTrials & TE.reject == 0 & ismember(TE.filename, midSessions);
lateTrials = csPlusTrials & bigRewardTrials & TE.reject == 0 & ismember(TE.filename, lateSessions);
%% make continuous early -> middle -> late phRaster

figSize = [1.66 0.68];
climfactor = 3;  
fdField = 'ZS';
linewidth = 4;
window = [-3 3];
offset = -0.1;
% tcolors = [204 153 255; 153 51 255; 51 0 102]; tcolors = tcolors ./ 255; 
tcolors = [0 0 0; 0.8 0.8 0.8; 0.6 0.6 0.6];
saveName = ['acquisition_phRaster_' animal];
ensureFigure(saveName, 1); axes; hold on;
phRasterFromTE(TE, csPlusTrials & bigRewardTrials & TE.reject == 0, 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry', 'zeroTimes', TE.usZeros, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,

% make lines to label early/middle/late trial sets for averages
allTrials = TE.filename(csPlusTrials & bigRewardTrials & TE.reject == 0);
theseTrials = find(ismember(allTrials, earlySessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(1,:), 'LineWidth', linewidth);
theseTrials = find(ismember(allTrials, midSessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(2,:), 'LineWidth', linewidth);
theseTrials = find(ismember(allTrials, lateSessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(3,:), 'LineWidth', linewidth);

% set(gca, 'Visible', 'off');
set(gca, 'YTick', [100 200], 'XTick', []);
colorbar;
formatFigurePublish('size', figSize, 'fontSize', 6);
if saveOn 
%     export_fig(fullfile(savepath, saveName), '-eps');
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));
end

%% make continuous early -> middle -> late lick raster

figSize = [1.66 0.68];


linewidth = 4;
window = [-3 3];
offset = -0.1;
% tcolors = [204 153 255; 153 51 255; 51 0 102]; tcolors = tcolors ./ 255; 
tcolors = [0 0 0; 0.8 0.8 0.8; 0.6 0.6 0.6];
saveName = ['acquisition_lickRaster_' animal];
ensureFigure(saveName, 1); axes; hold on;
eventRasterFromTE(TE, csPlusTrials & bigRewardTrials & TE.reject == 0, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroTimes', TE.usZeros, 'window', window); % 'CLimFactor', CLimFactor,

% make lines to label early/middle/late trial sets for averages
allTrials = TE.filename(csPlusTrials & bigRewardTrials & TE.reject == 0);
theseTrials = find(ismember(allTrials, earlySessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(1,:), 'LineWidth', linewidth);
theseTrials = find(ismember(allTrials, midSessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(2,:), 'LineWidth', linewidth);
theseTrials = find(ismember(allTrials, lateSessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(3,:), 'LineWidth', linewidth);

% set(gca, 'Visible', 'off');
set(gca, 'YTick', [100 200], 'XTick', [-3 0 3], 'XLim', window);
formatFigurePublish('size', figSize, 'fontSize', 6);

if saveOn 
%     export_fig(fullfile(savepath, saveName), '-eps');
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));
end

%% make overlaid phAverages

figSize = [1.66 0.9];
fdField = 'ZS';
% tcolors = [204 153 255; 153 51 255; 51 0 102]; tcolors = tcolors ./ 255; 
tcolors = [0 0 0; 0.8 0.8 0.8; 0.6 0.6 0.6];

saveName = ['overlaid_phAvg_' animal];
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, {earlyTrials, midTrials, lateTrials}, 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'alpha', 1, 'cmap', tcolors); % cued reward
set(gca, 'XLim', window, 'XTick', [-3 0 3], 'YTick', [0 2 4]);
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');
set(gca, 'XLim', window);
% xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize, 'fontSize', 6);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% averages by condition, steady state
window = [-3 3];
figSize = [1.48 0.76];
saveName = ['phAvg_complete_' animal];
ensureFigure(saveName, 1);
    trialSets = {...
        csPlusTrials & bigRewardTrials & ismember(TE.sessionIndex, 3:6);...
        uncuedTrials & bigRewardTrials & ismember(TE.sessionIndex, 3:6);...
        csPlusTrials & omitTrials & ismember(TE.sessionIndex, 3:6)...
        };
[ha, hl] = phPlotAverageFromTE(TE, trialSets, 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'linespec', {'b', 'c', 'k'}, 'alpha', 1); % cued reward
set(gca, 'XLim', window, 'XTick', [-3 0 3]);
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);
% xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize);
set(gca, 'YTick', [0 4]);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end    