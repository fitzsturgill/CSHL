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

%% signal correlations

xData = [us_pooled.rew.avg_mean(:,1) us_pooled.puff.avg_mean(:,1) us_pooled.shock.avg_mean(:,1)];
yData = [us_pooled.rew.avg_mean(:,2) us_pooled.puff.avg_mean(:,2) us_pooled.shock.avg_mean(:,2)];


nAnimals = size(xData, 1);
nPerm = 10000;
xShuff = repmat(xData, nPerm, 1);
yShuff = zeros(size(xShuff));
for counter = 0:nPerm - 1
    pix = randperm(size(yData, 1));
    yShuff(1 + counter*nAnimals:nAnimals + counter*nAnimals, :) = yData(pix, :);
end

rs = zeros(size(xData, 2),1);

for counter = 1:size(xData, 1)
    rs(counter) = corr(xData(counter,:)', yData(counter, :)');
end

rs_shuff = zeros(size(xShuff, 1),1);
for counter = 1:size(xShuff,1)
    rs_shuff(counter) = corr(xShuff(counter,:)', yShuff(counter, :)');
end

Rsignal = cum(rs);
Rsignal_shuff = cum(rs_shuff);
%%
saveName = 'cumHist_Rsignal';
ensureFigure(saveName, 1);
axes; hold on;
plot(Rsignal.sorted, Rsignal.index, '-k');
plot(Rsignal_shuff.sorted, Rsignal_shuff.index, 'Color', [0.8 0.8 0.8]);
xlabel('Rsignal');
formatFigurePublish('size', [1 1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end