%{
Script to generate simple RPE graphs and early/middle/late acquisition graphs for figure 1 from
ChAT_22, derived from McKnight_Poster_script
Fig1_GCaMP_simple.m
%}

%%
DB = dbLoadExperiment('SO_RewardPunish_odor');
photometryField = 'Photometry';
saveOn = 1;
climfactor = 2;
window = [-4 4];

animal = 'ChAT_22';
success = dbLoadAnimal(DB, animal);

% set savepath
savepath = fullfile(DB.path, ['figure1' filesep animal]);
ensureDirectory(savepath);

earlySessions = {'ChAT_22_SO_RewardPunish_odor_May11_2016_Session1.mat'};
midSessions = {'ChAT_22_SO_RewardPunish_odor_May12_2016_Session1.mat'};
lateSessions = {'ChAT_22_SO_RewardPunish_odor_May15_2016_Session1.mat', 'ChAT_22_SO_RewardPunish_odor_May16_2016_Session1.mat'};
lateLickSessions = {'ChAT_22_SO_RewardPunish_odor_May15_2016_Session1.mat'}; % May 16th session is "funny" for licking, i.e. lick sensor not working properly or mouse having incredibly well timed and efficient licks

%% make early/middle/late phRasters

figSize = [1.6 0.5];
climfactor = 3;  
fdField = 'ZS';

tcolors = [204 153 255; 153 51 255; 51 0 102]; tcolors = tcolors ./ 255; 

saveName = ['early_phRaster_' animal];
ensureFigure(saveName, 1);
phRasterFromTE(TE, rewardOdorTrials & rewardTrials & ismember(TE.filename, earlySessions), 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry', 'zeroTimes', TE.usZeros, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
% set(gca, 'Visible', 'off');
set(gca, 'YTick', [1 round(2/3 * sum(rewardOdorTrials & rewardTrials & ismember(TE.filename, earlySessions))/10)*10], 'XTick', []);
formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end    

saveName = ['middle_phRaster_' animal];
ensureFigure(saveName, 1);
phRasterFromTE(TE, rewardOdorTrials & rewardTrials & ismember(TE.filename, midSessions), 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry', 'zeroTimes', TE.usZeros, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
% set(gca, 'Visible', 'off');
set(gca, 'YTick', [1 round(2/3 * sum(rewardOdorTrials & rewardTrials & ismember(TE.filename, midSessions))/10)*10], 'XTick', []);
formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

saveName = ['late_phRaster_' animal];
ensureFigure(saveName, 1);
phRasterFromTE(TE, rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions), 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry', 'zeroTimes', TE.usZeros, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
% set(gca, 'Visible', 'off');
set(gca, 'YTick', [1 round(0.6 * sum(rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions))/10)*10], 'XTick', []);
formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end


%% make early/middle/late phAverages

figSize = [2 0.5];
fdField = 'ZS';
tcolor = mycolors('chat');

saveName = ['early_phAvg_' animal];
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, rewardOdorTrials & rewardTrials & ismember(TE.filename, earlySessions), 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'cmap', tcolor, 'alpha', 1); % cued reward
set(gca, 'YLim', [-1 2.5], 'XTickLabel', {},  'YTickLabel', {});
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end    

saveName = ['middle_phAvg_' animal];
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, rewardOdorTrials & rewardTrials & ismember(TE.filename, midSessions), 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'cmap', tcolor, 'alpha', 1); % cued reward
set(gca, 'YLim', [-1 2.5], 'XTickLabel', {},  'YTickLabel', {});
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end    

saveName = ['late_phAvg_' animal];
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions), 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'cmap', tcolor, 'alpha', 1); % cued reward
set(gca, 'YLim', [-1 2.5], 'XTickLabel', {},  'YTickLabel', {});
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end    

%% make overlaid phAverages

figSize = [1 1];
fdField = 'ZS';
tcolors = [204 153 255; 153 51 255; 51 0 102]; tcolors = tcolors ./ 255; 

saveName = ['overlaid_phAvg_' animal];
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, {rewardOdorTrials & rewardTrials & ismember(TE.filename, earlySessions),...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, midSessions),...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions)}, 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'alpha', 1, 'cmap', tcolors); % cued reward
set(gca, 'YLim', [-1 2.5],  'YTickLabel', {});
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);
xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize, 'fontSize', 6);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end    


%% phAverages, cued, uncued, omission

figSize = [2 1];
saveName = ['phAvg_complete_' animal];
ensureFigure(saveName, 1);
trialSets = {...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    uncuedTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    rewardOdorTrials & ~rewardTrials & ismember(TE.filename, lateSessions)...
    };
[ha, hl] = phPlotAverageFromTE(TE, trialSets, 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'linespec', {'b', 'c', 'k'}, 'alpha', 1); % cued reward
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);
xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end    

%% lick averages, cued, uncued, omission

figSize = [2 1];
tcolor = mycolors('chat');

saveName = ['lickAvgs_complete_' animal];
ensureFigure(saveName, 1);

varargin = {'trialNumbering', 'consecutive',...
    'window', [-1 4], 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'b', 'c', 'k'}, 'alpha', 1};

[ha, hl] = plotEventAverageFromTE(TE, {...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    uncuedTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    rewardOdorTrials & ~rewardTrials & ismember(TE.filename, lateSessions)}, 'Port1In', varargin{:});
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
ylabel('Licks (1/s)');  set(gca, 'XLim', window);
xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end    