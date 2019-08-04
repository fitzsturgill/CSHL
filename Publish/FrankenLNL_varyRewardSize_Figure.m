%{
FrankenLNL_varyRewardSize_Figure


%}

DB = dbLoadExperiment('FrankenLNL_varyRewardSize');
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);
saveOn = 1;

load(fullfile(DB.path, 'pooled', 'grandAverages.mat'));
disp(['*** loading: ' fullfile(DB.path, 'pooled', 'grandAverages.mat') ' ***']);

load(fullfile(DB.path, 'pooled', 'us_pooled.mat'));
disp(['*** loading: ' fullfile(DB.path, 'pooled', 'us_pooled.mat') ' ***']);


%% run stats


comp = {'large_vs_medium'; 'medium_vs_small'};
ph_n = zeros(2,1);
ph_p = zeros(2,1);
licks_p = zeros(2,1);
licks_n = zeros(2,1);
rewardSize_stats = table(comp, ph_p, ph_n, licks_p, licks_n);
rewardSize_stats.ph_p(1) = signrank(us_pooled.large.phPeakMean(:), us_pooled.medium.phPeakMean(:));
rewardSize_stats.ph_p(2) = signrank(us_pooled.medium.phPeakMean(:), us_pooled.small.phPeakMean(:));
rewardSize_stats.ph_n(:) = numel(us_pooled.large.phPeakMean);
rewardSize_stats.licks_p(1) = signrank(us_pooled.large.lickRate(:), us_pooled.medium.lickRate(:));
rewardSize_stats.licks_p(2) = signrank(us_pooled.medium.lickRate(:), us_pooled.small.lickRate(:));
rewardSize_stats.licks_n(:) = numel(us_pooled.large.lickRate);

if saveOn
    save(fullfile(savepath, 'rewardSize_stats.mat'), 'rewardSize_stats');    
end

rewardSize_stats

%% plot and bar graph

figSize = [1 1];
saveName = 'varyRewardSize_pooled_barGraph';
ensureFigure(saveName, 1);
data = [us_pooled.large.phPeakMean(:) us_pooled.medium.phPeakMean(:) us_pooled.small.phPeakMean(:)];
errorbar(3:-1:1, mean(data), std(data) / size(data, 1), 'k');
set(gca, 'YLim', [0.75 2]);
set(gca, 'XLim', [.5 3.5], 'XTick', [1 2 3], 'XTickLabel', fliplr({'10uL', '5uL', '2uL'}), 'YTick', [1 2]);
xlabel('Rew. size (uL)'); ylabel('Fluor. (\sigma-bl.)');

formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
end