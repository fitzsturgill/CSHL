% FrankenLNL_RewardPunish_correlations
% run FrankenLNL_RewardPunish_poolAnimals first

% plot noise correlations, etc. for FrankenLNL_RewardPunish


DB = dbLoadExperiment('FrankenLNL_RewardPunish');
saveOn = 1;

savepath = fullfile(DB.path, ['pooled' filesep]);
ensureDirectory(savepath);
figsavepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(figsavepath);


%% load data
load(fullfile(savepath, 'us_pooled.mat'), 'us_pooled');
disp(['*** loading: ' fullfile(savepath, 'us_pooled.mat') ' ***']);


%% plot cum hists of noise correlations expressed as Rsq (fraction of variance explained)
rewData = us_pooled.rew.Rnoise_mean;
puffData = us_pooled.puff.Rnoise_mean;
shockData = us_pooled.shock.Rnoise_mean;
blData = nanmean([us_pooled.rew.Rnoise_bl us_pooled.puff.Rnoise_bl us_pooled.shock.Rnoise_bl], 2);
Rnoise.reward = cum(rewData);
Rnoise.puff = cum(puffData);
Rnoise.shock = cum(shockData);
Rnoise.bl = cum(blData);
Rnoise.all = cum([rewData; puffData; shockData; blData]);

saveName = 'cumHist_Rnoise';
ensureFigure(saveName, 1);
axes; hold on;
plot(Rnoise.bl.sorted, Rnoise.bl.index, '-k');
plot(Rnoise.reward.sorted, Rnoise.reward.index, '-b');
plot(Rnoise.puff.sorted, Rnoise.puff.index, '-r');
plot(Rnoise.shock.sorted, Rnoise.shock.index, 'Color', mycolors('shock'));
xlabel('Rnoise');
% legend({'bl.', 'rew.', 'puff', 'shock'}, 'Location', 'northwest'); legend boxoff;
formatFigurePublish('size', [1 1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

saveName = 'cumHist_Rnoise_all';
ensureFigure(saveName, 1);
axes; hold on;
plot(Rnoise.all.sorted, Rnoise.all.index, '-k');
    