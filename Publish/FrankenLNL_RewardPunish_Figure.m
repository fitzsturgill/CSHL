%



DB = dbLoadExperiment('FrankenLNL_RewardPunish');
figPath = fullfile(DB.path, 'figure');
ensureDirectory(figPath);
smoothWindow = 1;
saveOn = 1;
figSize = [2 1];

%% examples for ACh_7 and ACh_15, showing similar signals, but different degrees of noise correlations

photometryField = 'Photometry';
fdField = 'ZS';
window = [-6.5 4];
linecolors = [0 0 1; 0 1 1];
lineWidth = 0.25;

animal = 'ACh_7';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups
saveName = sprintf('PE_BLA_%s_avgs', animal);  
h=ensureFigure(saveName, 1); 
sessionIndexList = 5; % just the 1 session because early sessions surprise modulation hasn't developed whereas later sessions I shorten the delay (surprise modulation is consistent) but it messes up the graph having 2 different delays


avgData1 = phAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 1,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window);
avgData2 = phAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 2,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window);
rch = [max(range(avgData1.Avg')) max(range(avgData2.Avg'))];
offset = rch(1) + rch(2) * 0.25;
xData = avgData1.xData(1,:);

subplot(1,2,2); hold on;
[thisHl, thisHp] = boundedline(xData, avgData1.Avg, permute(avgData1.SEM, [2 3 1]), gca, 'cmap', linecolors, 'alpha', 'nan', 'gap');
set(thisHl, 'LineWidth', lineWidth);
[thisHl, thisHp] = boundedline(xData, avgData2.Avg + offset, permute(avgData2.SEM, [2 3 1]), gca, 'cmap', linecolors, 'alpha', 'nan', 'gap');
set(thisHl, 'LineWidth', lineWidth);
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
set(gca, 'XLim', window);
xlabel('Time (s)');

animal = 'ACh_15';
sessionIndexList = [6 7];
success = dbLoadAnimal(DB, animal); % load TE and trial lookups
avgData1 = phAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 1,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window);
avgData2 = phAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 2,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window);
avgData1.Avg = avgData1.Avg * 3;
avgData1.SEM = avgData1.SEM * 3;
rch = [max(range(avgData1.Avg')) max(range(avgData2.Avg'))];
offset = rch(1) + rch(2) * 0.25;
xData = avgData1.xData(1,:);

subplot(1,2,1); hold on;
[thisHl, thisHp] = boundedline(xData, avgData1.Avg, permute(avgData1.SEM, [2 3 1]), gca, 'cmap', linecolors, 'alpha', 'nan', 'gap');
set(thisHl, 'LineWidth', lineWidth);
[thisHl, thisHp] = boundedline(xData, avgData2.Avg + offset, permute(avgData2.SEM, [2 3 1]), gca, 'cmap', linecolors, 'alpha', 'nan', 'gap');
set(thisHl, 'LineWidth', lineWidth);
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
set(gca, 'XLim', window);
xlabel('Time (s)');
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');

formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
    saveas(gcf, fullfile(figPath, [saveName '.fig']));
    saveas(gcf, fullfile(figPath, [saveName '.jpg'])); 
end