% cuedOutcome_pooled_figure

% goals: asummary data for cue and outcome responses


DB = dbLoadExperiment('cuedOutcome');

photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
ch = 1;

savePath = fullfile(DB.path, 'pooled');
ensureDirectory(savePath);
figPath = fullfile(DB.path, 'figure');
ensureDirectory(savePath);
nAnimals = length(DB.animals);

%%
% compile cue and outcome responses for high value, low value, and null,
% reward and punish


s2 = struct(...
    'data', [],...
    'avg', zeros(nAnimals, 1),...
    'SEM', zeros(nAnimals, 1)...
    );
cs_pooled = struct(...
    'high', s2,...
    'low', s2,...
    'uncued', s2...
    );
us_pooled = struct(...
    'reward', s2,... % uncued
    'punish', s2...  % uncued
    );


for counter = 1:nAnimals
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);        
    
    thisData = TE.phPeak_cs.data(highValueTrials);
    cs_pooled.high.data{counter} = thisData;
    cs_pooled.high.avg(counter) = nanmean(thisData);
    cs_pooled.high.SEM(counter) = nanSEM(thisData);
    
    thisData = TE.phPeak_cs.data(lowValueTrials);
    cs_pooled.low.data{counter} = thisData;
    cs_pooled.low.avg(counter) = nanmean(thisData);
    cs_pooled.low.SEM(counter) = nanSEM(thisData);
    
    thisData = TE.phPeak_us.data(uncuedTrials & rewardTrials);
    us_pooled.reward.data{counter} = thisData;
    us_pooled.reward.avg(counter) = nanmean(thisData);
    us_pooled.reward.SEM(counter) = nanSEM(thisData);
    
    thisData = TE.phPeak_us.data(uncuedTrials & rewardTrials);
    us_pooled.reward.data{counter} = thisData;
    us_pooled.reward.avg(counter) = nanmean(thisData);
    us_pooled.reward.SEM(counter) = nanSEM(thisData);    
    
    thisData = TE.phPeak_us.data(uncuedTrials & punishTrials);
    us_pooled.punish.data{counter} = thisData;
    us_pooled.punish.avg(counter) = nanmean(thisData);
    us_pooled.punish.SEM(counter) = nanSEM(thisData);        
end

%% also gather baselined us responses

s2 = struct(...
    'data', cell(1,1),...
    'avg', zeros(nAnimals, 1),...
    'SEM', zeros(nAnimals, 1)...
    );
cue = struct(...
    'high', s2,...
    'low', s2,...
    'cued', s2,...
    'uncued', s2...
    );
us_pooled_bl = struct(...
    'reward', cue,...
    'punish', cue,...
    'omit', cue...
    );


for counter = 1:nAnimals
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);        
        
    ordering = {...
        'reward', 'high', highValueTrials & rewardTrials;...
        'reward', 'low', lowValueTrials & rewardTrials;...
        'reward', 'cued', ~uncuedTrials & rewardTrials;...
        'reward', 'uncued', uncuedTrials & rewardTrials;...
        'punish', 'high', highValueTrials & punishTrials;...
        'punish', 'low', lowValueTrials & punishTrials;...
        'punish', 'cued', ~uncuedTrials & punishTrials;...        
        'punish', 'uncued', uncuedTrials & punishTrials;...        
        'omit', 'high', highValueTrials & omitTrials;...
        'omit', 'low', lowValueTrials & omitTrials;...
        'omit', 'cued', ~uncuedTrials & omitTrials;...        
        'omit', 'uncued', uncuedTrials & omitTrials;...                
        };
    
        for c2 = 1:size(ordering,1)
            thisData = TE.phPeak_us.data(ordering{c2, 3}) - TE.phPeak_preUs.data(ordering{c2, 3});
            us_pooled_bl.(ordering{c2,1}).(ordering{c2,2}).data{counter,:} = thisData;
            us_pooled_bl.(ordering{c2,1}).(ordering{c2,2}).avg(counter,:) = nanmean(thisData);
            us_pooled_bl.(ordering{c2,1}).(ordering{c2,2}).SEM(counter,:) = nanSEM(thisData);            
        end
    
end

%%
save(fullfile(savePath, 'us_pooled.mat'), 'us_pooled');
disp(['*** saving: ' fullfile(savePath, 'us_pooled.mat') ' ***']);
save(fullfile(savePath, 'cs_pooled.mat'), 'cs_pooled');
disp(['*** saving: ' fullfile(savePath, 'cs_pooled.mat') ' ***']);

%%
load(fullfile(savePath, 'us_pooled.mat'), 'us_pooled');
disp(['*** loading: ' fullfile(savePath, 'us_pooled.mat') ' ***']);
load(fullfile(savePath, 'cs_pooled.mat'), 'cs_pooled');
disp(['*** loading: ' fullfile(savePath, 'cs_pooled.mat') ' ***']);

%% us scatter plot

% Us scatter plot 
figSize = [1 1];
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\Adam\';
saveName = 'CuedOutcome_Us_scatterPlot';
ensureFigure(saveName, 1); 

% {0.9258, 0.4883, 0.1914} % orange
% {0.6680, 0.2148, 0.8359} % purple
% scatter(sumData.phReward_mean_bl.avg,sumData.phPunish_mean_bl.avg, 42, [0.6680, 0.2148, 0.8359], 'filled');
errorbar(us_pooled.reward.avg, us_pooled.punish.avg, -us_pooled.punish.SEM, us_pooled.punish.SEM,...
    -us_pooled.reward.SEM, us_pooled.reward.SEM, '.', 'Color', [0 0 0]);

% xlabel('Reward (\fontsize{20}\sigma\fontsize{16}-baseline)'); ylabel('Punish (\fontsize{20}\sigma\fontsize{16}-baseline)'); 
% xlabel('Reward (\sigma-baseline)'); ylabel('Punish (\sigma-baseline)'); 
ylim = get(gca, 'YLim'); xlim = get(gca, 'XLim');
set(gca, 'XTick', [0 1], 'YTick', [0 2], 'YLim', [-0.5 ylim(2)], 'XLim', [-0.5 xlim(2)]);
addOrginLines(gca, [0.7 0.7 0.7]);
xlabel('Reward'); ylabel('Punish');
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
    saveas(gcf, fullfile(figPath, [saveName '.fig']));
    saveas(gcf, fullfile(figPath, [saveName '.jpg']));   
end  
%%
[p, h] = ttest(us_pooled.reward.avg)

[p, h] = ttest(us_pooled.punish.avg)

%% cs scatter plot
    % cs scatter plot
saveName = 'CuedOutcome_Cs_scatterPlot';
ensureFigure(saveName, 1); 
figSize = [1 1.2];

scatter(cs_pooled.low.avg,cs_pooled.high.avg, 36, [0 0 0], 'filled');
errorbar(cs_pooled.low.avg,cs_pooled.high.avg, -cs_pooled.high.SEM, cs_pooled.high.SEM,...
    -cs_pooled.low.SEM, cs_pooled.low.SEM, '.', 'Color', [0 0 0]);
% xlabel('Low Value (\fontsize{20}\sigma\fontsize{16}-baseline)'); ylabel('High Value (\fontsize{20}\sigma\fontsize{16}-baseline)');

% set(gca, 'YLim', [0 1]); set(gca, 'XLim', [0 1]);
setXYsameLimit;
set(gca, 'XTick', [0 1], 'YTick', [0 1]);
addUnityLine(gca, [0.7 0.7 0.7]);
xlabel('Low value'); ylabel('High value');
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
    saveas(gcf, fullfile(figPath, [saveName '.fig']));
    saveas(gcf, fullfile(figPath, [saveName '.jpg']));   
end  

[p, h] = ttest(cs_pooled.low.avg, cs_pooled.high.avg)
%% us bl scatter plot

% Us scatter plot 
figSize = [1 1.2];
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\Adam\';
saveName = 'delta_Us_scatterPlot';
ensureFigure(saveName, 1); 

% scatter(sumData.phReward_mean_bl.avg,sumData.phPunish_mean_bl.avg, 42, [0.6680, 0.2148, 0.8359], 'filled');
errorbar(us_pooled_bl.reward.low.avg, us_pooled_bl.reward.high.avg, -us_pooled_bl.reward.high.SEM, us_pooled_bl.reward.high.SEM,...
    -us_pooled_bl.reward.low.SEM, us_pooled_bl.reward.low.SEM, '.', 'Color', [0 0 0]);

% xlabel('Reward (\fontsize{20}\sigma\fontsize{16}-baseline)'); ylabel('Punish (\fontsize{20}\sigma\fontsize{16}-baseline)'); 
% xlabel('Reward (\sigma-baseline)'); ylabel('Punish (\sigma-baseline)'); 
% ylim = get(gca, 'YLim'); xlim = get(gca, 'XLim');
% set(gca, 'XTick', [0 1], 'YTick', [0 2], 'YLim', [-0.5 ylim(2)], 'XLim', [-0.5 xlim(2)]);
% addOrginLines(gca, [0.7 0.7 0.7]);
setXYsameLimit;
addUnityLine(gca, [0.7 0.7 0.7]);
set(gca, 'XTick', [0 1], 'YTick', [0 1]);
xlabel('Low value'); ylabel('High value');
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
    saveas(gcf, fullfile(figPath, [saveName '.fig']));
    saveas(gcf, fullfile(figPath, [saveName '.jpg']));   
end  


%% us bl bar plot
% Us scatter plot 
figSize = [1.2 1.1];
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\Adam\';
saveName = 'delta_Us_barPlot';
ensureFigure(saveName, 1);
yData = [nanmean(us_pooled_bl.reward.uncued.avg) nanmean(us_pooled_bl.reward.low.avg) nanmean(us_pooled_bl.reward.high.avg)];
bData = [nanSEM(us_pooled_bl.reward.uncued.avg) nanSEM(us_pooled_bl.reward.low.avg) nanSEM(us_pooled_bl.reward.high.avg)];

errorbar([1 2 3], yData, bData, 'k');
set(gca, 'XLim', [0.5 3.5], 'YLim', [0 1.1], 'XTickLabel', {'', 'low', 'high'});

formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end

%% grand averages aligned to cue and licking

s2 = struct(...
    'xData', [],...
    'Avg', [],...
    'SEM', []...
    );
cue = struct(...
    'high', s2,...
    'low', s2,...
    'cued', s2,...
    'uncued', s2...
    );
outcome = struct(...
    'reward', cue,... % uncued
    'punish', cue,...  % uncued
    'omit', cue...
    );

[cueAligned, lickAligned, outcomes] = deal(outcome);


w1 = [-3 3]; %old: [-2 5]
w2 = [-3 3]; %old: [-2 5]
wPre = [-0.5 0];
wOutcome = [-0.5 3];

for counter = 1:nAnimals
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);
    
    ordering = {...
        'reward', 'high', highValueTrials & rewardTrials;...
        'reward', 'low', lowValueTrials & rewardTrials;...
        'reward', 'cued', ~uncuedTrials & rewardTrials;...
        'reward', 'uncued', uncuedTrials & rewardTrials;...
        'punish', 'high', highValueTrials & punishTrials;...
        'punish', 'low', lowValueTrials & punishTrials;...
        'punish', 'cued', ~uncuedTrials & punishTrials;...        
        'punish', 'uncued', uncuedTrials & punishTrials;...        
        'omit', 'high', highValueTrials & omitTrials;...
        'omit', 'low', lowValueTrials & omitTrials;...
        'omit', 'cued', ~uncuedTrials & omitTrials;...        
        'omit', 'uncued', uncuedTrials & omitTrials;...                
        };
    
    for c2 = 1:size(ordering,1)
        v1 = {'FluorDataField', 'ZS', 'window', w1, 'zeroTimes', TE.Cue, 'PhotometryField', 'Photometry'};
        avgData = phAverageFromTE(TE, ordering{c2,3}, 1, v1{:});
        cueAligned.(ordering{c2,1}).(ordering{c2,2}).Avg(counter,:) = avgData.Avg;
        cueAligned.(ordering{c2,1}).(ordering{c2,2}).SEM(counter,:) = avgData.SEM;
        cueAligned.(ordering{c2,1}).(ordering{c2,2}).xData(counter,:) = avgData.xData;
        
        v2 = {'FluorDataField', 'ZS', 'window', w2, 'zeroTimes', TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue), 'PhotometryField', 'Photometry'};
        avgData = phAverageFromTE(TE, ordering{c2,3}, 1, v2{:});
        lickAligned.(ordering{c2,1}).(ordering{c2,2}).Avg(counter,:) = avgData.Avg;
        lickAligned.(ordering{c2,1}).(ordering{c2,2}).SEM(counter,:) = avgData.SEM;
        lickAligned.(ordering{c2,1}).(ordering{c2,2}).xData(counter,:) = avgData.xData;
        
        vPre = {'FluorDataField', 'ZS', 'window', wPre, 'zeroTimes', TE.Us, 'PhotometryField', 'Photometry'};
        avgData_pre = phAverageFromTE(TE, ordering{c2,3}, 1, vPre{:});        
        
        vOutcome = {'FluorDataField', 'ZS', 'window', wOutcome, 'zeroTimes', TE.Us, 'PhotometryField', 'Photometry'};
        avgData_outcome = phAverageFromTE(TE, ordering{c2,3}, 1, vOutcome{:});
        outcomes.(ordering{c2,1}).(ordering{c2,2}).Avg(counter,:) = avgData_outcome.Avg - nanmean(avgData_pre.Avg);
        outcomes.(ordering{c2,1}).(ordering{c2,2}).SEM(counter,:) = avgData_outcome.SEM;
        outcomes.(ordering{c2,1}).(ordering{c2,2}).xData(counter,:) = avgData_outcome.xData;        
    end
end

%% example rasters cue and lick aligned to accompany averages
animal = 'ChAT_42';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups

figSize = [1.1 1.5];
photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
%
channels = [1];
climfactor = 2;
lickTickLineWidth = 0.3; % 0.3 works well when printed out given current sizing
markerSize = 1;
fontsize = 10;


xwindow = [-3 3];


lickOnsets = TE.lickLatency_cs(highValueTrials & rewardTrials);
lickOnsets = sort(lickOnsets);    
lickZeros = TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue);

saveName = 'lickAligned_highValue_rasters_left';  
ensureFigure(saveName, 1); 
eventRasterFromTE(TE, highValueTrials & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.lickLatency_cs, 'LineWidth', lickTickLineWidth);
set(gca, 'XLim', xwindow);
set(gca, 'YTickLabel', {});
set(gca, 'XTick', [-3 0 3]);
xlabel('Time from odor (s)');
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

saveName = 'lickAligned_highValue_rasters_middle';  
ensureFigure(saveName, 1); 
phRasterFromTE(TE, highValueTrials & rewardTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'zeroTimes', TE.Cue, 'window', xwindow, 'sortValues', TE.lickLatency_cs); % 'CLimFactor', CLimFactor,
line(lickOnsets, (1:sum(highValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');    
set(gca, 'YTickLabel', {});
set(gca, 'XTick', [-3 0 3], 'XLim', xwindow);
xlabel('Time from odor (s)');
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

saveName = 'lickAligned_highValue_rasters_right';  
ensureFigure(saveName, 1); 
axes; hold on;
phRasterFromTE(TE, highValueTrials & rewardTrials, 1, 'zeroTimes', lickZeros, 'window', xwindow,...
    'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'showSessionBreaks', 0, 'sortValues', TE.lickLatency_cs); % 'CLimFactor', CLimFactor,
line(-1 * lickOnsets, (1:sum(highValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');    
% scatter(-1 * TE.lickLatency_cs(highValueTrials & rewardTrials), 1:sum(highValueTrials & rewardTrials), markerSize, [1 1 1], 'filled');
xlabel('Time from lick (s)');
set(gca, 'YTickLabel', {});
set(gca, 'XTick', [-3 0 3], 'XLim', xwindow);
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  



lickOnsets = TE.lickLatency_cs(lowValueTrials & rewardTrials);
lickOnsets = sort(lickOnsets);    
lickZeros = TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue);

saveName = 'lickAligned_lowValue_rasters_left';  
ensureFigure(saveName, 1); 
eventRasterFromTE(TE, lowValueTrials & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.lickLatency_cs, 'LineWidth', lickTickLineWidth);
set(gca, 'XLim', xwindow);
set(gca, 'YTickLabel', {});
xlabel('Time from odor (s)');
set(gca, 'XTick', [-3 0 3]);
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

saveName = 'lickAligned_lowValue_rasters_middle';  
ensureFigure(saveName, 1); 
phRasterFromTE(TE, lowValueTrials & rewardTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'zeroTimes', TE.Cue, 'window', xwindow, 'sortValues', TE.lickLatency_cs); % 'CLimFactor', CLimFactor,
line(lickOnsets, (1:sum(lowValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');    
set(gca, 'YTickLabel', {});
xlabel('Time from odor (s)');
set(gca, 'XTick', [-3 0 3], 'XLim', xwindow);
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

saveName = 'lickAligned_lowValue_rasters_right';  
ensureFigure(saveName, 1); 
axes; hold on;
phRasterFromTE(TE, lowValueTrials & rewardTrials, 1, 'zeroTimes', lickZeros, 'window', xwindow,...
    'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'showSessionBreaks', 0, 'sortValues', TE.lickLatency_cs); % 'CLimFactor', CLimFactor,
line(-1 * lickOnsets, (1:sum(lowValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');    
xlabel('Time from lick (s)');
set(gca, 'YTickLabel', {});
set(gca, 'XTick', [-3 0 3], 'XLim', xwindow);
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

%% corresponding averages

ylim = [-1 2];
figSize = [1.66 0.5];
saveName = 'lickAligned_highValue_avgs_middle';  
ensureFigure(saveName, 1); 
phPlotAverageFromTE(TE, highValueTrials & rewardTrials, 1, 'zeroTimes', TE.Cue, 'window', xwindow, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'cmap', [1 0 1]);
set(gca, 'XTick', [0], 'YTick', [0 1], 'YLim', ylim, 'YTickLabel', {}, 'XTickLabel', {});
addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

saveName = 'lickAligned_highValue_avgs_right';  
ensureFigure(saveName, 1); 
phPlotAverageFromTE(TE, highValueTrials & rewardTrials, 1, 'zeroTimes', TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue), 'window', xwindow, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'cmap', [1 0 1]);
set(gca, 'XTick', [0], 'YTick', [0 1], 'YLim', ylim, 'YTickLabel', {}, 'XTickLabel', {});
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

saveName = 'lickAligned_lowValue_avgs_middle';  
ensureFigure(saveName, 1); 
phPlotAverageFromTE(TE, lowValueTrials & rewardTrials, 1, 'zeroTimes', TE.Cue, 'window', xwindow, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'cmap', [0 0 1]);
set(gca, 'XTick', [0], 'YTick', [0 1], 'YLim', ylim, 'YTickLabel', {}, 'XTickLabel', {});
addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

saveName = 'lickAligned_lowValue_avgs_right';  
ensureFigure(saveName, 1); 
phPlotAverageFromTE(TE, lowValueTrials & rewardTrials, 1, 'zeroTimes', TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue), 'window', xwindow, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'cmap', [0 0 1]);
set(gca, 'XTick', [0], 'YTick', [0 1], 'YLim', ylim, 'YTickLabel', {}, 'XTickLabel', {});
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end

%% just high value, reward aligned, for figure 1 to match time scale of spike data

animal = 'ChAT_42';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups

figSize = [1.7 1];
photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
%
channels = [1];
climfactor = 2;
lickTickLineWidth = 0.3; % 0.3 works well when printed out given current sizing
markerSize = 1;
fontsize = 10;


xwindow = [-4 3];


lickOnsets = TE.lickLatency_cs(highValueTrials & rewardTrials) - 3; % rezero to reward
lickOnsets = sort(lickOnsets);    
% lickZeros = TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue);

saveName = 'highValue_phRaster_lickSorted_figure1';  
ensureFigure(saveName, 1); 
phRasterFromTE(TE, highValueTrials & rewardTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'zeroTimes', TE.Us, 'window', xwindow, 'sortValues', TE.lickLatency_cs); % 'CLimFactor', CLimFactor,
line(lickOnsets, (1:sum(highValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-');    
set(gca, 'YTickLabel', {}, 'XTick', [-3 0 3], 'YTick', [0 200 400], 'XLim', xwindow);
% xlabel('Time from odor (s)');
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

ylim = [-1 2];
figSize = [1.66 0.5];a
saveName = 'highValue_phAverages_lickSorted_figure1';  
ensureFigure(saveName, 1); 
phPlotAverageFromTE(TE, highValueTrials & rewardTrials, 1, 'zeroTimes', TE.Us, 'window', xwindow, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'cmap', [1 0 1]);
set(gca, 'XTick', [0], 'YTick', [0 1], 'YLim', ylim, 'YTickLabel', {}, 'XTickLabel', {});
addStimulusPatch(gca, [-3 -2], '', [0.8 0.8 0.8], 0.5);
addStimulusPatch(gca, [-0.1 0.1], '', [0.8 0.8 0.8], 0.5);
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
end  

%% subtract off cue-aligned mean just for fun

saveName = 'test';  
ensureFigure(saveName, 1); 
axes; hold on;

ts = highValueTrials & rewardTrials;
predata = TE.Photometry.data.ZS(ts, :);
predata = predata - nanmean(predata);
cdata = alignedDataWindow(predata, true(size(predata, 1), 1), 'zeroTimes', lickZeros(ts), 'window', xwindow, 'Fs', 20, 'startTimes', TE.Photometry.startTime(ts));
imagesc(cdata, [-2 2]); % 'CLimFactor', CLimFactor,

xlabel('Time from lick (s)');
set(gca, 'YTickLabel', {});
% formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
end  


%%


lickOnsets = TE.lickLatency_cs(lowValueTrials & rewardTrials);
lickOnsets = sort(lickOnsets);        
subplot(1,4,3);
eventRasterFromTE(TE, lowValueTrials & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.lickLatency_cs, 'LineWidth', lickTickLineWidth);
set(gca, 'XLim', xwindow);
% title('low value');

subplot(1,4,4);
phRasterFromTE(TE, lowValueTrials & rewardTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'sortValues', TE.lickLatency_cs, 'zeroTimes', TE.Cue, 'window', xwindow); % 'CLimFactor', CLimFactor,             
line(lickOnsets, (1:sum(lowValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1);
axs = findobj(gcf, 'Type', 'axes');
set(axs, 'FontSize', fontsize);
set(axs([1 3]), 'YTick', []);

formatFigurePublish('size', [4 1.2]);

if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
%     export_fig(fullfile(savePath, saveName), '-eps');
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));   
end  



%% plot grand averages conditioned on cue responses

saveName = 'grandAverages_cuedOutcome';
ensureFigure(saveName, 1);
figSize = [1.1 1.5];
ylim = [-0.2 2];
axes;
yData = [nanmean(cueAligned.reward.high.Avg)' nanmean(cueAligned.reward.low.Avg)'];
bData = [nanSEM2(cueAligned.reward.high.Avg)' nanSEM2(cueAligned.reward.low.Avg)'];
bData = permute(bData, [1 3 2]);
cmap = [1 0 1; 0 0 1];
boundedline(cueAligned.reward.cued.xData(1,:)', yData, bData, 'cmap', cmap);
addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);
addStimulusPatch(gca, [2.9 3.1], '', [0.8 0.8 0.8], 0.5); 
% xlabel('time from odor (s)'); ylabel('Fluor. (\sigma-bl.)');
set(gca, 'XLim', w1, 'YLim', ylim, 'XTick', [-3 0 3]);
formatFigurePublish('size', figSize);    

if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
end

saveName = 'grandAverages_cuedOutcome_lickAligned';
ensureFigure(saveName, 1);

axes;
yData = [nanmean(lickAligned.reward.high.Avg)' nanmean(lickAligned.reward.low.Avg)'];
bData = [nanSEM2(lickAligned.reward.high.Avg)' nanSEM2(lickAligned.reward.low.Avg)'];
bData = permute(bData, [1 3 2]);
cmap = [1 0 1; 0 0 1];
boundedline(lickAligned.reward.cued.xData(1,:)', yData, bData, 'cmap', cmap);
addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);
addStimulusPatch(gca, [2.9 3.1], '', [0.8 0.8 0.8], 0.5); 
% xlabel('time from first lick (s)'); ylabel('Fluor. (\sigma-bl.)');
set(gca, 'XLim', w1, 'YLim', ylim, 'XTick', [-3 0 3]);
formatFigurePublish('size', figSize);    

if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
end

%% re-baselined outcome responses to argue for surprise modulation from GCaMP


saveName = 'grandAverages_reward_bl';
ensureFigure(saveName, 1);
figSize = [2 1];
% ylim = [-0.2 2];
axes;
yData = [nanmean(outcomes.reward.high.Avg)' nanmean(outcomes.reward.low.Avg)' nanmean(outcomes.reward.uncued.Avg)'];
bData = [nanSEM2(outcomes.reward.high.Avg)' nanSEM2(outcomes.reward.low.Avg)' nanSEM2(outcomes.reward.uncued.Avg)'];
bData = permute(bData, [1 3 2]);
cmap = [1 0 1; 0 0 1; 0 0 0];
boundedline(outcomes.reward.cued.xData(1,:)', yData, bData, 'cmap', cmap);
addStimulusPatch(gca, [-0.1 0.1], '', [0.8 0.8 0.8], 0.5);
xlabel('Time from reward (s)'); ylabel('Fluor. (\sigma-bl.)');
% set(gca, 'XLim', w1, 'YLim', ylim);
formatFigurePublish('size', figSize);    

if saveOn
%     print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
end

saveName = 'grandAverages_punish_bl';
ensureFigure(saveName, 1);
figSize = [2 1];
% ylim = [-0.2 2];
axes;
yData = [nanmean(outcomes.punish.high.Avg)' nanmean(outcomes.punish.low.Avg)' nanmean(outcomes.punish.uncued.Avg)'];
bData = [nanSEM2(outcomes.punish.high.Avg)' nanSEM2(outcomes.punish.low.Avg)' nanSEM2(outcomes.punish.uncued.Avg)'];
bData = permute(bData, [1 3 2]);
cmap = [1 0 1; 0 0 1; 0 0 0];
boundedline(outcomes.punish.cued.xData(1,:)', yData, bData, 'cmap', cmap);
addStimulusPatch(gca, [-0.1 0.1], '', [0.8 0.8 0.8], 0.5);
xlabel('Time from punish(s)'); ylabel('Fluor. (\sigma-bl.)');
% set(gca, 'XLim', w1, 'YLim', ylim);
formatFigurePublish('size', figSize);    

if saveOn
%     print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
end



%% Hacky deconvolved version:
% deconvolved

% goals: summary data for cue and outcome responses


DB = dbLoadExperiment('cuedOutcome');

photometryField = 'Photometry';
fdField = 'ZSdeconv';
saveOn = 1;
ch = 1;

savePath = fullfile(DB.path, 'pooled');
ensureDirectory(savePath);
figPath = fullfile(DB.path, 'figure');
ensureDirectory(savePath);
nAnimals = length(DB.animals);


% grand averages aligned to cue and licking

s2 = struct(...
    'xData', [],...
    'Avg', [],...
    'SEM', []...
    );
cue = struct(...
    'high', s2,...
    'low', s2,...
    'cued', s2,...
    'uncued', s2...
    );
outcome = struct(...
    'reward', cue,... % uncued
    'punish', cue,...  % uncued
    'omit', cue...
    );

[cueAligned, lickAligned] = deal(outcome);


w1 = [-2 5];
w2 = [-2 5];

for counter = 1:nAnimals
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);
    
    ordering = {...
        'reward', 'high', highValueTrials & rewardTrials;...
        'reward', 'low', lowValueTrials & rewardTrials;...
        'reward', 'cued', ~uncuedTrials & rewardTrials;...
        'reward', 'uncued', uncuedTrials & rewardTrials;...
        'punish', 'high', highValueTrials & punishTrials;...
        'punish', 'low', lowValueTrials & punishTrials;...
        'punish', 'cued', ~uncuedTrials & punishTrials;...        
        'punish', 'uncued', uncuedTrials & punishTrials;...        
        'omit', 'high', highValueTrials & omitTrials;...
        'omit', 'low', lowValueTrials & omitTrials;...
        'omit', 'cued', ~uncuedTrials & omitTrials;...        
        'omit', 'uncued', uncuedTrials & omitTrials;...                
        };
    
    for c2 = 1:size(ordering,1)
        v1 = {'FluorDataField', fdField, 'window', w1, 'zeroTimes', TE.Cue, 'PhotometryField', 'Photometry'};
        avgData = phAverageFromTE(TE, ordering{c2,3}, 1, v1{:});
        cueAligned.(ordering{c2,1}).(ordering{c2,2}).Avg(counter,:) = avgData.Avg;
        cueAligned.(ordering{c2,1}).(ordering{c2,2}).SEM(counter,:) = avgData.SEM;
        cueAligned.(ordering{c2,1}).(ordering{c2,2}).xData(counter,:) = avgData.xData;
        
        v2 = {'FluorDataField', fdField, 'window', w2, 'zeroTimes', TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue), 'PhotometryField', 'Photometry'};
        avgData = phAverageFromTE(TE, ordering{c2,3}, 1, v2{:});
        lickAligned.(ordering{c2,1}).(ordering{c2,2}).Avg(counter,:) = avgData.Avg;
        lickAligned.(ordering{c2,1}).(ordering{c2,2}).SEM(counter,:) = avgData.SEM;
        lickAligned.(ordering{c2,1}).(ordering{c2,2}).xData(counter,:) = avgData.xData;
    end
end



% plot grand averages conditioned on cue responses

saveName = 'grandAverages_cuedOutcome_deconv';
ensureFigure(saveName, 1);
figSize = [2 1];
ylim = [-0.2 2];
axes;
yData = [nanmean(cueAligned.reward.high.Avg)' nanmean(cueAligned.reward.low.Avg)'];
bData = [nanSEM2(cueAligned.reward.high.Avg)' nanSEM2(cueAligned.reward.low.Avg)'];
bData = permute(bData, [1 3 2]);
cmap = [1 0 1; 0 0 1];
boundedline(cueAligned.reward.cued.xData(1,:)', yData, bData, 'cmap', cmap);
addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);
addStimulusPatch(gca, [2.9 3.1], '', [0.8 0.8 0.8], 0.5); 
xlabel('time from odor (s)'); ylabel('Fluor. (\sigma-bl.)');
set(gca, 'XLim', w1, 'YLim', ylim);
formatFigurePublish('size', figSize);    

if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
end

saveName = 'grandAverages_cuedOutcome_lickAligned_deconv';
ensureFigure(saveName, 1);

axes;
yData = [nanmean(lickAligned.reward.high.Avg)' nanmean(lickAligned.reward.low.Avg)'];
bData = [nanSEM2(lickAligned.reward.high.Avg)' nanSEM2(lickAligned.reward.low.Avg)'];
bData = permute(bData, [1 3 2]);
cmap = [1 0 1; 0 0 1];
boundedline(lickAligned.reward.cued.xData(1,:)', yData, bData, 'cmap', cmap);
addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);
addStimulusPatch(gca, [2.9 3.1], '', [0.8 0.8 0.8], 0.5); 
xlabel('time from first lick (s)'); ylabel('Fluor. (\sigma-bl.)');
set(gca, 'XLim', w1, 'YLim', ylim);
formatFigurePublish('size', figSize);    

if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
end