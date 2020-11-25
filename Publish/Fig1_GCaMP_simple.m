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
window = [-3 3];

animal = 'ChAT_22';
success = dbLoadAnimal(DB, animal);

% set savepath
savepath = fullfile(DB.path, ['figure1' filesep animal]);
ensureDirectory(savepath);

earlySessions = {'ChAT_22_SO_RewardPunish_odor_May11_2016_Session1.mat'};
midSessions = {'ChAT_22_SO_RewardPunish_odor_May12_2016_Session1.mat'};
lateSessions = {'ChAT_22_SO_RewardPunish_odor_May15_2016_Session1.mat', 'ChAT_22_SO_RewardPunish_odor_May16_2016_Session1.mat'};
lateLickSessions = {'ChAT_22_SO_RewardPunish_odor_May15_2016_Session1.mat'}; % May 16th session is "funny" for licking, i.e. lick sensor not working properly or mouse having incredibly well timed and efficient licks
fdField = 'ZS';
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
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
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
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
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
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
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
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
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
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
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
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end    

%% make overlaid phAverages

figSize = [1.66 0.9];
fdField = 'ZS';
% tcolors = [204 153 255; 153 51 255; 51 0 102]; tcolors = tcolors ./ 255; 
tcolors = [0 0 0; 0.8 0.8 0.8; 0.6 0.6 0.6];

saveName = ['overlaid_phAvg_' animal];
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, {rewardOdorTrials & rewardTrials & ismember(TE.filename, earlySessions),...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, midSessions),...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions)}, 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'alpha', 1, 'cmap', tcolors); % cued reward
set(gca, 'YLim', [-1 2.5], 'XLim', window,  'YTickLabel', {}, 'XTick', [-3 0 3]);
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);
% xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize, 'fontSize', 6);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% make overlaid lick averages

figSize = [1.66 0.9];

% tcolors = [204 153 255; 153 51 255; 51 0 102]; tcolors = tcolors ./ 255; 
tcolors = [0 0 0; 0.8 0.8 0.8; 0.6 0.6 0.6];

saveName = ['overlaid_lickAvg_' animal];
ensureFigure(saveName, 1);
[ha, hl] = plotEventAverageFromTE(TE, {rewardOdorTrials & rewardTrials & ismember(TE.filename, earlySessions),...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, midSessions),...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, lateLickSessions)}, 'Port1In',...
    'zeroTimes', TE.usZeros, 'window', window, 'alpha', 1, 'cmap', tcolors); % cued reward
set(gca, 'XLim', window, 'XTick', [-3 0 3]);
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);
% xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize, 'fontSize', 6);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%     export_fig(fullfile(savepath, saveName), '-eps');
end

%% phRasters, cued, uncued, omission, example #1
figSize = [1.6 0.81];
window = [-3 3];
tcolor = mycolors('chat');
trialSets = {...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    uncuedTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    rewardOdorTrials & ~rewardTrials & ismember(TE.filename, lateSessions)...
    };

saveName= ['phRasters_complete_' animal '_cued'];
ensureFigure(saveName, 1);
climfactor = 3;  
axes;
phRasterFromTE(TE, rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions), 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', 'ZS', 'PhotometryField', 'Photometry', 'zeroTimes', TE.Reward, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
set(gca, 'YTick', [1 50 100], 'YTickLabel', {'1', '50', ''});
formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end        

saveName= ['phRasters_complete_' animal '_uncued'];
ensureFigure(saveName, 1);
climfactor = 3;  
axes;
phRasterFromTE(TE, uncuedTrials & rewardTrials & ismember(TE.filename, lateSessions), 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', 'ZS', 'PhotometryField', 'Photometry', 'zeroTimes', TE.Reward, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
set(gca, 'YTick', [1 50 100], 'YTickLabel', {'1', '50', ''});
formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end       


%% phAverages, cued, uncued, omission, example #1

figSize = [1.48 0.76];
saveName = ['phAvg_complete_' animal];
ensureFigure(saveName, 1);
trialSets = {...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    uncuedTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    rewardOdorTrials & ~rewardTrials & ismember(TE.filename, lateSessions)...
    };
[ha, hl] = phPlotAverageFromTE(TE, trialSets, 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'linespec', {'b', 'c', 'k'}, 'alpha', 1); % cued reward
set(gca, 'XLim', window, 'XTick', [-3 0 3]);
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);
% xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end    

%% lick averages, cued, uncued, omission

figSize = [1.48 0.76];
tcolor = mycolors('chat');

saveName = ['lickAvgs_complete_' animal];
ensureFigure(saveName, 1);

varargin = {'trialNumbering', 'consecutive',...
    'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'b', 'c', 'k'}, 'alpha', 1};

[ha, hl] = plotEventAverageFromTE(TE, {...
    rewardOdorTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    uncuedTrials & rewardTrials & ismember(TE.filename, lateSessions),...
    rewardOdorTrials & ~rewardTrials & ismember(TE.filename, lateSessions)}, 'Port1In', varargin{:});
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
set(gca, 'XLim', window, 'XTick', [-3 0 3], 'YTick', [0 10]);
% xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%     export_fig(fullfile(savepath, saveName), '-eps');
end    

%% phAverages, cued, uncued, omission, example #2
animal = 'ChAT_26';
success = dbLoadAnimal(DB, animal);

% set savepath
savepath = fullfile(DB.path, ['figure1' filesep animal]);
ensureDirectory(savepath);

figSize = [1.48 0.76];
saveName = ['phAvg_complete_' animal];
ensureFigure(saveName, 1);
trialSets = {...
    rewardOdorTrials & rewardTrials,...
    uncuedTrials & rewardTrials,...
    rewardOdorTrials & ~rewardTrials...
    };
[ha, hl] = phPlotAverageFromTE(TE, trialSets, 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'window', window, 'linespec', {'b', 'c', 'k'}, 'alpha', 1); % cued reward
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);
xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end    

%% aversive averages
animal = 'ChAT_20';
success = dbLoadAnimal(DB, animal);
window = [-4 4];
% set savepath
savepath = fullfile(DB.path, ['figure1' filesep animal]);
ensureDirectory(savepath);


figSize = [1.48 0.76];
saveName = ['phAvg_aversive_' animal];
ensureFigure(saveName, 1);
trialSets = {...
    punishOdorTrials & punishTrials,...
    uncuedTrials & punishTrials,...
    punishOdorTrials & ~punishTrials...
    };
[ha, hl] = phPlotAverageFromTE(TE, trialSets, 1,...
    'zeroTimes', TE.usZeros, 'FluorDataField', fdField, 'linespec', {'r', 'm', 'k'}, 'window', window, 'alpha', 1); % cued reward , 'linespec', {'r', 'm', 'k'}
set(gca, 'XLim', window, 'XTick', [-3 0 3]);
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);
% xlabel('Time from reinforcement (s)');
formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%     export_fig(fullfile(savepath, saveName), '-eps');
end  