% FrankenLNL_RewardPunish_correlations
% run FrankenLNL_RewardPunish_poolAnimals first

% plot noise correlations, etc. for FrankenLNL_RewardPunish


DB = dbLoadExperiment('FrankenLNL_RewardPunish');
saveOn = 1;

savePath = fullfile(DB.path, ['pooled' filesep]);
ensureDirectory(savePath);
figPath = fullfile(DB.path, 'figure');
ensureDirectory(figPath);


%% load data
load(fullfile(savePath, 'us_pooled.mat'), 'us_pooled');
disp(['*** loading: ' fullfile(savePath, 'us_pooled.mat') ' ***']);

%% load data
load(fullfile(savePath, 'cs_pooled.mat'), 'cs_pooled');
disp(['*** loading: ' fullfile(savePath, 'cs_pooled.mat') ' ***']);


%% plot cum hists of noise correlations expressed as R

% point out 2 example mice with vertical dotted lines
figSize = [2 1];
examples = {'ACh_15', 'ACh_7'};
ix = find(ismember(DB.animals, examples));

rewData = us_pooled.rew.Rnoise_mean;
puffData = us_pooled.puff.Rnoise_mean;
shockData = us_pooled.shock.Rnoise_mean;
rewData_cued = us_pooled.rew_cued.Rnoise_mean;
puffData_cued = us_pooled.puff_cued.Rnoise_mean;
shockData_cued = us_pooled.shock_cued.Rnoise_mean;
csPlusData = cs_pooled.CSplus.Rnoise_mean;
csMinusData = cs_pooled.CSminus.Rnoise_mean;

exData_rew = repmat(rewData(ix)', 2, 1);
exData_rewCued = repmat(rewData_cued(ix)', 2, 1);
exData_y = repmat([0; 1], 1, size(exData_rew, 2));

blData = nanmean([us_pooled.rew.Rnoise_bl us_pooled.puff.Rnoise_bl us_pooled.shock.Rnoise_bl us_pooled.rew_cued.Rnoise_bl us_pooled.puff_cued.Rnoise_bl us_pooled.shock_cued.Rnoise_bl], 2);
Rnoise.reward = cum(rewData);
Rnoise.puff = cum(puffData);
Rnoise.shock = cum(shockData);
Rnoise.reward_cued = cum(rewData_cued);
Rnoise.puff_cued = cum(puffData_cued);
Rnoise.shock_cued = cum(shockData_cued);
Rnoise.bl = cum(blData);
Rnoise.all = cum([rewData; puffData; shockData; rewData_cued; puffData_cued; shockData_cued; blData]);
Rnoise.csPlus = cum(csPlusData);
Rnoise.csMinus = cum(csMinusData);

saveName = 'cumHist_Rnoise';
ensureFigure(saveName, 1);
axes; hold on;
line(exData_rew, exData_y, 'Color', [0 1 1], 'LineStyle', '--');
line(exData_rewCued, exData_y, 'Color', [0 0 1], 'LineStyle', '--');
plot(Rnoise.bl.sorted, Rnoise.bl.index, 'Color', 'k');
plot(Rnoise.reward.sorted, Rnoise.reward.index, 'Color', 'c');
plot(Rnoise.puff.sorted, Rnoise.puff.index, 'Color', 'm');
plot(Rnoise.shock.sorted, Rnoise.shock.index, 'Color', 'g');
plot(Rnoise.reward_cued.sorted, Rnoise.reward_cued.index, 'Color', 'b');
plot(Rnoise.puff_cued.sorted, Rnoise.puff_cued.index, 'Color', 'r');
plot(Rnoise.shock_cued.sorted, Rnoise.shock_cued.index, 'Color', [0 0.5 0]);
plot(Rnoise.csPlus.sorted, Rnoise.csPlus.index, ':b');
plot(Rnoise.csMinus.sorted, Rnoise.csMinus.index, ':r');
xlabel('Rnoise');
set(gca, 'XLim', [-1 1], 'YTick', [0 1]);
% legend({'bl.', 'rew.', 'puff', 'shock'}, 'Location', 'northwest'); legend boxoff;
formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
    saveas(gcf, fullfile(figPath, [saveName '.fig']));
    saveas(gcf, fullfile(figPath, [saveName '.jpg']));
end

saveName = 'cumHist_Rnoise_all';
ensureFigure(saveName, 1);
axes; hold on;
plot(Rnoise.all.sorted, Rnoise.all.index, '-k');

%% signal correlations

xData = [us_pooled.rew.avg_mean(:,1) us_pooled.puff.avg_mean(:,1) us_pooled.shock.avg_mean(:,1) us_pooled.rew_cued.avg_mean(:,1) us_pooled.puff_cued.avg_mean(:,1) us_pooled.shock_cued.avg_mean(:,1) cs_pooled.CSplus.avg_mean(:,1) cs_pooled.CSminus.avg_mean(:,1)];
yData = [us_pooled.rew.avg_mean(:,2) us_pooled.puff.avg_mean(:,2) us_pooled.shock.avg_mean(:,2) us_pooled.rew_cued.avg_mean(:,2) us_pooled.puff_cued.avg_mean(:,2) us_pooled.shock_cued.avg_mean(:,2) cs_pooled.CSplus.avg_mean(:,2) cs_pooled.CSminus.avg_mean(:,2)];


nAnimals = size(xData, 1);
nPerm = 100;
xShuff = repmat(xData, nPerm, 1);
yShuff = zeros(size(xShuff));
for counter = 0:nPerm - 1
    pix = randperm(size(yData, 1));
    yShuff(1 + counter*nAnimals:nAnimals + counter*nAnimals, :) = yData(pix, :);
end

[rs, rsq, slopes1, slopes2] = deal(zeros(size(xData, 2),1));

% line fit
ft = fittype({'x'});
for counter = 1:size(xData, 1)
    fitx = xData(counter,:)';
    fity = yData(counter,:)';
    rs(counter) = corr(fitx, fity);
    [fo, gof, output] = fit(fitx, fity, ft);
    slopes1(counter) = fo.a;
    cvd = cov(fitx, fity);
    cvd = cvd(2); % off diagonal
    slope = covData + mean(fitx) * mean(fity) / (var(fitx) + mean(fitx)^2);
    slopes2(counter) = slope;
    rsq(counter) = gof.rsquare;
end

[rs_shuff, rsq_shuff] = deal(zeros(size(xShuff, 1),1));
for counter = 1:size(xShuff,1)
    fitx = xShuff(counter,:)';
    fity = yShuff(counter, :)';
    rs_shuff(counter) = corr(fitx, fity);
    [fo, gof, output] = fit(fitx, fity, ft);
    rsq_shuff(counter) = gof.rsquare;
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

Rsq_signal = cum(rsq);
Rsq_signal_shuff = cum(rsq_shuff);
saveName = 'cumHist_Rsq_LineThroughOrigin_signal';
ensureFigure(saveName, 1);
axes; hold on;
plot(Rsq_signal.sorted, Rsq_signal.index, '-k');
plot(Rsq_signal_shuff.sorted, Rsq_signal_shuff.index, 'Color', [0.8 0.8 0.8]);
xlabel('Rsq');
formatFigurePublish('size', [1 1]);

if saveOn 
    export_fig(fullfile(savePath, saveName), '-eps');
end
%% show distributions and means in an array for all mice

FluorField_us = 'phPeakMean_us';

saveName = sprintf('LeftRight_Distributions_%s', FluorField_us);
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
    linecolors = [0 1 1; 1 0 1; 0 1 0; 0 0 1; 1 0 0; 0 0.5 0];         
    trialSets = {'rew', 'puff', 'shock', 'rew_cued', 'puff_cued', 'shock_cued'};
    trialSetNames = {'Reward', 'Air Puff', 'Shock', 'Reward cued', 'Puff cued', 'Shock cued'};   
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
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
end


%% show example distribution

FluorField_us = 'phPeakMean_us';

saveName = sprintf('LeftRight_Distributions_example_%s', FluorField_us);
ensureFigure(saveName, 1);
mcPortraitFigSetup(gcf);
% formatFigurePublish('size', [4 2], 'fontSize', 8);
% animals = {'ACh_7', 'ACh_3'};
animals = {'ACh_7'};



% xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];

nCol = ceil(sqrt(length(animals)));
nRow = ceil(length(animals) / nCol);
linecolors = [0 1 1; 1 0 1; 0 1 0; 0 0 1; 1 0 0; 0 0.5 0];         
trialSets = {'rew', 'puff', 'shock', 'rew_cued', 'puff_cued', 'shock_cued'};
trialSetNames = {'Reward', 'Air Puff', 'Shock', 'Reward cued', 'Puff cued', 'Shock cued'};   
for acounter = 1:length(animals)
    subplot(nCol, nRow,acounter); hold on;        
    % reward is first in the trialSets list, use it to normalize

    xDataStruct = struct('Mean', cell(size(trialSets, 2), 1), 'SEM', cell(size(trialSets, 2), 1));
    yDataStruct = struct('Mean', cell(size(trialSets, 2), 1), 'SEM', cell(size(trialSets, 2), 1));
    h=[];   
    [xMean, yMean] = deal(zeros(length(trialSets), 1));
    for counter = 1:length(trialSets)       
        setField = trialSets{counter};
        xData = us_pooled.(setField).mean{acounter}(:,1); yData = us_pooled.(setField).mean{acounter}(:,2); 
        h(end + 1) = scatter(xData, yData, 2, linecolors(counter, :), '.', 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 0.2);%, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
        xMean(counter) = nanmean(xData);
        yMean(counter) = nanmean(yData);
    end

    h = [];
    for counter = 1:size(trialSets, 2)
        plot([0 xMean(counter)], [0 yMean(counter)], 'Color', linecolors(counter, :), 'Linewidth', 2);
        
        scatter(xMean(counter), yMean(counter), 40, linecolors(counter, :), 'o', 'MarkerFaceAlpha', 0, 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', 'none');            
%         h(end + 1) = plot(xMean(counter), yMean(counter), '-o', 'Color', linecolors(counter, :));
    end
    

%     sameXYScale(gca);
    addOrginLines(gca);
%     legend(h, trialSetNames, 'Location', 'Best'); legend('boxoff'); 
    % subplot(1,2,1);
    % textBox('Cue', [], [0.5 0.95], 8);
    % xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
%     ylabel('\fontsize{8}Right (\fontsize{12}\sigma\fontsize{8}-baseline)');
    ylabel('Right');
    % subplot(1,2,2);
    % textBox('Outcome', [], [0.5 0.95], 8);
%     xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
    xlabel('Left');

    % ylabel('');   
end


formatFigurePublish('size', [1.5 1.5]);

if saveOn 
    export_fig(fullfile(savePath, saveName), '-eps');
end



%% calculate auROC
% assuming that a linear relationship describes response of both L and R
% BLA ACh release to different strengths of the different reinforcers/signals, then
% use L/R ratios as estimates the slopes of those linear relationships.
alpha = 0.01;
nBoot = 10000;
nAnimals = length(us_pooled.animals);
slopes = struct(...
    'auROC', zeros(nAnimals, 3),...
    'p', zeros(nAnimals, 3)...
    );
slopes.comparisons = {...
    'rew_vs_puff', 'rew', 'puff';...
    'rew_vs_shock', 'rew', 'shock';...
    'puff_vs_shock', 'puff', 'shock';...
    };

for acounter = 1:nAnimals
    for ccounter = 1:length(slopes.comparisons)
%         xData = us_pooled.(slopes.comparisons{ccounter, 2}).mean{acounter}(:,2) ./ us_pooled.(slopes.comparisons{ccounter, 2}).mean{acounter}(:,1); % rise over run
%         yData = us_pooled.(slopes.comparisons{ccounter, 3}).mean{acounter}(:,2) ./ us_pooled.(slopes.comparisons{ccounter, 3}).mean{acounter}(:,1); % rise over run
        
        
        % mean centered, idiot check, still results in non-zero auROC
        % values and rejected null hypotheses, what gives
        xData = (us_pooled.(slopes.comparisons{ccounter, 2}).mean{acounter}(:,2) - mean(us_pooled.(slopes.comparisons{ccounter, 2}).mean{acounter}(:,2))) ./ (us_pooled.(slopes.comparisons{ccounter, 2}).mean{acounter}(:,1) - mean(us_pooled.(slopes.comparisons{ccounter, 2}).mean{acounter}(:,1))); % rise over run
        yData = (us_pooled.(slopes.comparisons{ccounter, 3}).mean{acounter}(:,2) - mean(us_pooled.(slopes.comparisons{ccounter, 3}).mean{acounter}(:,2))) ./ (us_pooled.(slopes.comparisons{ccounter, 3}).mean{acounter}(:,1) - mean(us_pooled.(slopes.comparisons{ccounter, 3}).mean{acounter}(:,1))); % rise over run
        
        
        [D, P] = rocarea(xData, yData, 'boot', nBoot, 'scale');
        slopes.auROC(acounter, ccounter) = D;
        slopes.p(acounter, ccounter) = P;
    end
end
%
saveName = 'slopes_auROC_histogram';
ensureFigure(saveName, 1);
subplot(1,2,1); histogram(slopes.auROC(:), -1:0.1:1); xlabel('auROC scaled for slopes');
subplot(1,2,2); histogram(slopes.p(:), 0:0.01:1); xlabel('pvals for same');

binEdges = -1:0.2:1;
binCenters = binEdges(1:end-1) + (binEdges(2) - binEdges(1));
saveName = 'slopes_auROC_histogram2';
ensureFigure(saveName, 1);
axes; hold on;
% let's do the Benjamini-Hochberg stepdown procedure to control false
% discovery rate to be < alpha

% rejected = slopes.auROC(slopes.p < 0.05);
% accepted = slopes.auROC(slopes.p >= 0.05);
auROC_all = slopes.auROC(:);
p_all = slopes.p(:);
[p_sorted, I] = sort(p_all);
auROC_sorted = auROC_all(I);
alpha_adjusted = (1:nAnimals * 3)' ./ (nAnimals * 3) * alpha;
h_sorted = p_sorted < alpha_adjusted;
[rejected_counts, ~] = histcounts(auROC_sorted(h_sorted), binEdges);
[accepted_counts, ~] = histcounts(auROC_sorted(~h_sorted), binEdges);

b = bar(binCenters, [accepted_counts'  rejected_counts'], 'stacked');
b(1).FaceColor = 'w';
b(2).FaceColor = 'k';


%% lets fit the mean centered us responses for each mouse to a line, then use the slope to align the us means across reinforcements.  visualization of how divergent are the us responses

% fit all responses together with a linefit through the origin and use this
% for scaling
nAnimals = length(DB.animals);



ensureFigure('show_noise_fits', 1);
mcPortraitFigSetup(gcf);
nCol = ceil(sqrt(nAnimals));
nRow = ceil(nAnimals / nCol);

all_slopes = zeros(nAnimals, 1);

trialSets = {'rew', 'puff', 'shock', 'rew_cued', 'puff_cued', 'shock_cued'};
for counter = 1:nAnimals
    [ch1Data, ch2Data] = deal([]);
    for tscounter = 1:length(trialSets)
        ts = trialSets{tscounter};
        ch1Data = [ch1Data; us_pooled.(ts).mean{counter}(:,1)];% - us_pooled.(ts).avg_mean(counter,1)];
        ch2Data = [ch2Data; us_pooled.(ts).mean{counter}(:,2)];% - us_pooled.(ts).avg_mean(counter,2)];
    end
    % line fit
    model = 'a*x'; % has to pass through origin
    ft = fittype(model);
    fo = fitoptions('Method', 'LinearLeastSquares',...
        'Upper', Inf,...
        'Lower', -Inf ...
    );
    fo = fit(ch1Data, ch2Data, ft);
    all_slopes(counter) = fo.a;
    subplot(nCol, nRow, counter); hold on;
    scatter(ch1Data, ch2Data, 20, 'k', '.');
    plot(fo, 'predfunc'); legend off;
    setXYsameLimit;
    addOrginLines;
end

%

linecolors = [0 1 1; 1 0 1; 0 1 0; 0 0 1; 1 0 0; 0 0.5 0];         
trialSetNames = {'Reward', 'Air Puff', 'Shock', 'Reward cued', 'Puff cued', 'Shock cued'};   



[xData, yData] = deal(zeros(nAnimals, length(trialSets)));

for counter = 1:nAnimals
    for tscounter = 1:length(trialSets)
        ts = trialSets{tscounter};
        xData(counter, tscounter) = us_pooled.(ts).avg_mean(counter, 2) / all_slopes(counter);
        yData(counter, tscounter) = us_pooled.(ts).avg_mean(counter, 2);
    end
end

saveName = 'us_aligned_all';
ensureFigure(saveName, 1);
axes; hold on;
for tscounter = 1:length(trialSets)
    disp(num2str(tscounter))
    scatter(xData(:,tscounter), yData(:,tscounter), 20, linecolors(tscounter, :), 'o', 'filled');
end
addUnityLine;    


%% show distributions and means in an array for all mice

FluorField_us = 'phPeakMean_us';

saveName = sprintf('LeftRight_Distributions_%s', FluorField_us);
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
    linecolors = [0 1 1; 1 0 1; 0 1 0; 0 0 1; 1 0 0; 0 0.5 0];         
    trialSets = {'rew', 'puff', 'shock', 'rew_cued', 'puff_cued', 'shock_cued'};
    trialSetNames = {'Reward', 'Air Puff', 'Shock', 'Reward cued', 'Puff cued', 'Shock cued'};   
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
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
end

%% show distributions and means in an array for all mice- VERSION 2 with cue

FluorField_us = 'phPeakMean_us';

saveName = sprintf('LeftRight_DistributionsV2_%s', FluorField_us);
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
    linecolors = [0 0 1; 1 0 0; 0 0.5 0; 0 1 1; 1 0 1];         
    trialSets = {'rew_cued', 'puff_cued', 'shock_cued', 'CSplus', 'CSminus'};
    trialSetNames = {'Reward cued', 'Puff cued', 'Shock cued', 'cs+', 'cs-'};   
    xDataStruct = struct('Mean', cell(size(trialSets, 2), 1), 'SEM', cell(size(trialSets, 2), 1));
    yDataStruct = struct('Mean', cell(size(trialSets, 2), 1), 'SEM', cell(size(trialSets, 2), 1));
    h=[];   
    for counter = 1:length(trialSets)        
    %     allTrials = allTrials | trialSets{counter};
        setField = trialSets{counter};
        if ismember(setField, {'CSplus', 'CSminus'})
            xData = cs_pooled.(setField).mean{acounter}(:,1); yData = cs_pooled.(setField).mean{acounter}(:,2);
        else
            xData = us_pooled.(setField).mean{acounter}(:,1); yData = us_pooled.(setField).mean{acounter}(:,2); 
        end
        
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
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
end

    


