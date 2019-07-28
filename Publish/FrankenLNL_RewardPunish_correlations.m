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
%
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
%% 
%% Development code below: normalize responses to punishment by those to reward and compare between left and right BLA


FluorField_us = 'phPeakMean_us';

saveName = sprintf('LeftRight_PunNorm_%s', FluorField_us);
ensureFigure(saveName, 1);
mcPortraitFigSetup(gcf);
% formatFigurePublish('size', [4 2], 'fontSize', 8);
% animals = {'ACh_7', 'ACh_3'};
animals = DB.animals;


% xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];

nCol = ceil(sqrt(length(animals)));
nRow = ceil(length(animals) / nCol);

for acounter = 1:length(animals)
    subplot(nCol, nRow,acounter); hold on;    
    title(animals{acounter}, 'Interpreter', 'none');
    % reward is first in the trialSets list, use it to normalize
    linecolors = [0 0 1; 1 0 0; 0 1 0];         
    trialSets = {'rew', 'puff', 'shock'};
    trialSetNames = {'Reward', 'Air Puff', 'Shock'};   
    xDataStruct = struct('Mean', cell(size(trialSets, 2), 1), 'SEM', cell(size(trialSets, 2), 1));
    yDataStruct = struct('Mean', cell(size(trialSets, 2), 1), 'SEM', cell(size(trialSets, 2), 1));
    h=[];   
    for counter = 1:length(trialSets)        
    %     allTrials = allTrials | trialSets{counter};
        setField = trialSets{counter};
        xData = us_pooled.(setField).mean{acounter}(:,1); yData = us_pooled.(setField).mean{acounter}(:,2); 
        
%         if counter == 1
%             xDenom = nanmean(xData);
%             yDenom = nanmean(yData);
%         end
%         xData = xData ./ xDenom;
%         yData = yData ./ yDenom;

        h(end + 1) = scatter(xData, yData, 20, linecolors(counter, :), '.', 'MarkerFaceColor', 'flat');%, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);

        xDataStruct(counter).Mean  = nanmean(xData);
        xDataStruct(counter).SEM = nanstd(xData) / sqrt(sum(isfinite(xData)));
        yDataStruct(counter).Mean = nanmean(yData);
        yDataStruct(counter).SEM = nanstd(yData) / sqrt(sum(isfinite(yData)));
    end
%     addUnityLine(gca, [0.3 0.3 0.3]);

%     xlim = get(gca,'XLim');
%     ylim = get(gca,'YLim');
% 
%     p1 = min([min(xlim) min(ylim)]);
%     p2 = min([max(xlim) max(ylim)]);
%     p1 = 

%     h = line('Parent',gca,'XData',[p1 p2],'YData',[p2 p1]);

%     set(h,'Color',[0.3 0.3 0.3]);
    h = [];
    for counter = 1:size(trialSets, 2)
        plot([0 xDataStruct(counter).Mean], [0 yDataStruct(counter).Mean], 'Color', linecolors(counter, :), 'Linewidth', 3);
        
        errorbar(xDataStruct(counter).Mean, yDataStruct(counter).Mean, -yDataStruct(counter).SEM, yDataStruct(counter).SEM,-xDataStruct(counter).SEM, xDataStruct(counter).SEM, 'Color', linecolors(counter, :), 'LineWidth', 3, 'CapSize', 2, 'Marker', 'o');
        h(end + 1) = plot(xDataStruct(counter).Mean, yDataStruct(counter).Mean, '-', 'Color', linecolors(counter, :));
    end
    

%     sameXYScale(gca);
    addOrginLines(gca);
    legend(h, trialSetNames, 'Location', 'Best'); legend('boxoff'); 
    % subplot(1,2,1);
    % textBox('Cue', [], [0.5 0.95], 8);
    % xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
%     ylabel('\fontsize{8}Right (\fontsize{12}\sigma\fontsize{8}-baseline)');
    ylabel('Right (reward-normalized)');
    % subplot(1,2,2);
    % textBox('Outcome', [], [0.5 0.95], 8);
%     xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
    xlabel('Left (reward-normalized)');

    % ylabel('');   
end





if saveOn 
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
end


%% calculate auROC