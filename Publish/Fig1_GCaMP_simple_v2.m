%{
Script to generate simple RPE graphs and early/middle/late acquisition graphs for figure 1 from
ChAT_22, derived from McKnight_Poster_script
Fig1_GCaMP_simple_v2.m
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
phRasterFromTE(TE, rewardOdorTrials & rewardTrials, 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry', 'zeroTimes', TE.usZeros, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,

% make lines to label early/middle/late trial sets for averages
allTrials = TE.filename(rewardOdorTrials & rewardTrials);
theseTrials = find(ismember(allTrials, earlySessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(1,:), 'LineWidth', linewidth);
theseTrials = find(ismember(allTrials, midSessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(2,:), 'LineWidth', linewidth);
theseTrials = find(ismember(allTrials, lateSessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(3,:), 'LineWidth', linewidth);

% set(gca, 'Visible', 'off');
set(gca, 'YTick', [100 300], 'XTick', []);
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
eventRasterFromTE(TE, rewardOdorTrials & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroTimes', TE.usZeros, 'window', window); % 'CLimFactor', CLimFactor,

% make lines to label early/middle/late trial sets for averages
allTrials = TE.filename(rewardOdorTrials & rewardTrials);
theseTrials = find(ismember(allTrials, earlySessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(1,:), 'LineWidth', linewidth);
theseTrials = find(ismember(allTrials, midSessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(2,:), 'LineWidth', linewidth);
theseTrials = find(ismember(allTrials, lateSessions));
line([window(2) + offset window(2) + offset], [theseTrials(1) theseTrials(end)], 'Color', tcolors(3,:), 'LineWidth', linewidth);

% set(gca, 'Visible', 'off');
set(gca, 'YTick', [100 300], 'XTick', [-3 0 3]);
formatFigurePublish('size', figSize, 'fontSize', 6);
if saveOn 
%     export_fig(fullfile(savepath, saveName), '-eps');
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%     saveas(gcf, fullfile(savepath, [saveName '.epsc']));
end


%% align raster/ averages to odor onset vs lick onset, version #1

TE.lickLatency_cs = calcEventLatency(TE, 'Port1In', TE.Cue, TE.usZeros);
lickZeros = TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue);
specialWindowByLick = [cellfun(@(x) x(1), TE.Cue) cellfun(@(x) x(end), TE.Cue)] - lickZeros;

% rasters
saveName = ['align_lick_vs_odor_v1_' animal];
ensureFigure(saveName, 1); 


subplot(2,2,1);
specialWindowByCue = [0 2];
phRasterFromTE(TE, rewardOdorTrials, 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry',...
    'zeroTimes', TE.Cue, 'window', specialWindowByCue, 'sortValues', TE.lickLatency_cs, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
ylabel ('trial # sorted by lick latency');

subplot(2,2,2);
phRasterFromTE(TE, rewardOdorTrials, 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry',...
    'zeroTimes', lickZeros, 'window', specialWindowByLick, 'sortValues', TE.lickLatency_cs, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
 set(gca, 'YTick', []);

% averages
tcolor = mycolors('chat');
subplot(2,2,3); hold on;
phPlotAverageFromTE(TE, rewardOdorTrials, 1, 'FluorDataField', fdField, 'PhotometryField', 'Photometry',...
    'zeroTimes', TE.Cue, 'window', specialWindowByCue, 'cmap', tcolor, 'alpha', 1); % 'CLimFactor', CLimFactor,
xlabel('time from cue (s)'); 
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');

subplot(2,2,4);
phPlotAverageFromTE(TE, rewardOdorTrials, 1, 'FluorDataField', fdField, 'PhotometryField', 'Photometry',...
    'zeroTimes', lickZeros, 'window', specialWindowByLick, 'cmap', tcolor, 'alpha', 1); % 'CLimFactor', CLimFactor,
xlabel('time from response (s)');


%% align raster/ averages to odor onset vs lick onset, version #2

figSize = [2.5 1.7] * 1.2;

TE.lickLatency_cs = calcEventLatency(TE, 'Port1In', TE.Cue, TE.usZeros);
lickZeros = TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue);
specialWindow = [-1 1];

% rasters
saveName = ['align_lick_vs_odor_v2_' animal];
ensureFigure(saveName, 1); 


subplot(2,2,1); hold on;
specialWindowByCue = [0 2];
phRasterFromTE(TE, rewardOdorTrials, 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry',...
    'zeroTimes', TE.Cue, 'window', specialWindow, 'sortValues', TE.lickLatency_cs, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
lickTimes_sorted = sort(TE.lickLatency_cs(rewardOdorTrials));
plot(lickTimes_sorted, 1:sum(rewardOdorTrials), '--r');
ylabel ('trial # sorted'); set(gca, 'YTick', [100 300 500], 'XTick', []);

subplot(2,2,2); hold on;
phRasterFromTE(TE, rewardOdorTrials, 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry',...
    'zeroTimes', lickZeros, 'window', specialWindow, 'showSessionBreaks', 0, 'sortValues', TE.lickLatency_cs); % 'CLimFactor', CLimFactor,
plot(0 - lickTimes_sorted, 1:sum(rewardOdorTrials), '--k');
 set(gca, 'YTick', [], 'XTick', []);

% averages
tcolor = mycolors('chat');
ha=[];
ha(1) = subplot(2,2,3); hold on;
phPlotAverageFromTE(TE, rewardOdorTrials, 1, 'FluorDataField', fdField, 'PhotometryField', 'Photometry',...
    'zeroTimes', TE.Cue, 'window', specialWindow, 'cmap', tcolor, 'alpha', 1); % 'CLimFactor', CLimFactor,
xlabel('time from cue (s)'); 
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');

ha(2) = subplot(2,2,4); hold on;
phPlotAverageFromTE(TE, rewardOdorTrials, 1, 'FluorDataField', fdField, 'PhotometryField', 'Photometry',...
    'zeroTimes', lickZeros, 'window', specialWindow, 'cmap', tcolor, 'alpha', 1); % 'CLimFactor', CLimFactor,
xlabel('time from response (s)'); 

sameYScale(ha);
yData = get(gca, 'YLim');
plot([0 0], yData, '--r');
subplot(2,2,3);
plot([0 0], yData, '--k');

formatFigurePublish('size', figSize);
if saveOn 
%     export_fig(fullfile(savepath, saveName), '-eps');
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end
