%%
plotFields = {'blLicks', 'anticipatoryLicks1', 'anticipatoryLicks2', 'usLicks', 'epoch', 'sessionIndex',...
    'phCuePeak', 'phCuePeakLong', 'phDelayPeak', 'phUsPeak', 'trialIndex', 'phRaster'};
 
trialCutoffs = {...
    'ChAT_22_SO_RewardPunish_odor_May11_2016_Session1.mat', 200;...    
    'ChAT_22_SO_RewardPunish_odor_May12_2016_Session1.mat', 0;...
    'ChAT_22_SO_RewardPunish_odor_May13_2016_Session1.mat', 270;...
    'ChAT_22_SO_RewardPunish_odor_May14_2016_Session1.mat', 350;...
    'ChAT_22_SO_RewardPunish_odor_May15_2016_Session1.mat', 0;...
    };
includedTrials = zeros(size(trialCutoffs, 1), max(allTE.trialIndex));
for counter = 1:size(trialCutoffs, 1)
    if trialCutoffs{counter, 2} 
        includedTrials(counter, :) = bpFilterTrials2(allTE, 'fileName', trialCutoffs{counter, 1}, 'Trial', 1:trialCutoffs{counter, 2});
    else
        includedTrials(counter, :) = bpFilterTrials2(allTE, 'fileName', trialCutoffs{counter, 1});
    end
end
includedTrials = sum(includedTrials)'; % 'any' of each row of includedTrials
includedTrials(find(sum(isnan(allTE.phRaster), 2))) = 0; % get rid of NaN- containing trials


TrialTypes = [1:7];
byTypes = struct();
for counter = 1:length(TrialTypes)
    byTypes(counter).Trials = bpFilterTrials2(allTE, 'TrialTypes', TrialTypes(counter)) & includedTrials;
    for counter2 = 1:length(plotFields)
        thisData = allTE.(plotFields{counter2});
        if size(thisData, 1) == 1
            thisData = thisData';
        end
        byTypes(counter).(plotFields{counter2}) = thisData(byTypes(counter).Trials, :);
    end
end

%%
saveName = 'ChAT_22';
CLim = [-0.025 0.025];
bounds = [-3; 6];
types = [1 2 5];
xdata = linspace(bounds(1), bounds(2), size(allTE.phRaster, 2));
ensureFigure([saveName '_reward_phRasters'], 1);
ah = zeros(size(types));
for i = 1:length(types)
    type = types(i);
    ah(i) = subplot(1,length(types),i);
    imagesc(xdata, 1:size(byTypes(type).phRaster, 1), byTypes(type).phRaster); 
    breaks = find(diff(byTypes(type).sessionIndex))';
    line(repmat(bounds, 1, length(breaks)), [breaks; breaks], 'Parent', ah(i), 'Color', 'r'); % session breaks
end
set(ah, 'CLim', CLim, 'TickDir', 'out');

%%
% CLim = [-0.01 0.01];
% ensureFigure([saveName '_punish_phRasters'], 1); 
% ah = [];
% ah(end+1) = subplot(1,3,1);
% imagesc(xdata, 1:size(byTypes(3).phRaster, 1), byTypes(3).phRaster); 
% ah(end+1) = subplot(1,3,2);
% imagesc(xdata, 1:size(byTypes(4).phRaster, 1), byTypes(4).phRaster); 
% ah(end+1) = subplot(1,3,3);
% imagesc(xdata, 1:size(byTypes(6).phRaster, 1), byTypes(6).phRaster); 
% set(ah, 'CLim', CLim, 'TickDir', 'out')

%%
ensureFigure([saveName '_longitudinalCueResponses'], 1); 
subplot(4,1,1); scatter(byTypes(1).trialIndex, byTypes(1).anticipatoryLicks2 - byTypes(1).blLicks,'.'); % both blLicks and anticipatoryLicks2 are from 2s periods
% subplot(4,1,1); scatter(1:length(allTE.epoch), cuedReward.anticipatoryLicks2);
maxLR = max(byTypes(1).anticipatoryLicks2);
% subplot(4,1,1); plot(smooth(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
% maxLR = smooth(max(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
% hold on; stem(cuedReward.trialIndex(1:end-1), abs(diff(cuedReward.epoch)) * maxLR, 'g', 'Marker', 'none');
hold on; stem(byTypes(1).trialIndex(1:end-1), diff(byTypes(1).sessionIndex) * maxLR, 'r', 'Marker', 'none');
ylabel('delta lick rate (licks/s)'); title('cued reward trials');
set(gca, 'XLim', [1 length(allTE.epoch)]);

subplot(4,1,2); scatter(byTypes(1).trialIndex, byTypes(1).phCuePeak, '.');
maxP = max(byTypes(1).phCuePeak);
% hold on; stem(byTypes(1).trialIndex(1:end-1), abs(diff(cuedReward.epoch)) * maxP, 'g', 'Marker', 'none');
hold on; stem(byTypes(1).trialIndex(1:end-1), diff(byTypes(1).sessionIndex) * maxP, 'r', 'Marker', 'none');
ylabel('cue dF/F');
set(gca, 'XLim', [1 length(allTE.epoch)]); set(gca, 'YLim', [-0.02 0.02]);

subplot(4,1,3); scatter(byTypes(1).trialIndex, byTypes(1).phUsPeak, '.');
maxP = max(byTypes(1).phUsPeak);
% hold on; stem(byTypes(1).trialIndex(1:end-1), abs(diff(cuedReward.epoch)) * maxP, 'g', 'Marker', 'none');
hold on; stem(byTypes(1).trialIndex(1:end-1), diff(byTypes(1).sessionIndex) * maxP, 'r', 'Marker', 'none');
ylabel('cued US dF/F');
set(gca, 'XLim', [1 length(allTE.epoch)]); set(gca, 'YLim', [-0.02 0.02]);

subplot(4,1,4); scatter(byTypes(5).trialIndex, byTypes(5).phUsPeak, '.');
maxP = max(byTypes(5).phUsPeak);
% hold on; stem(byTypes(1).trialIndex(1:end-1), abs(diff(cuedReward.epoch)) * maxP, 'g', 'Marker', 'none');
hold on; stem(byTypes(5).trialIndex(1:end-1), diff(byTypes(5).sessionIndex) * maxP, 'r', 'Marker', 'none');
ylabel('uncued US dF/F');
set(gca, 'XLim', [1 length(allTE.epoch)]); set(gca, 'YLim', [-0.02 0.02]);

%%

cd('C:\Users\Adam\Dropbox\KepecsLab\_Fitz\McKnight\Poster\');
early_dFF_cued = byTypes(1).phRaster(byTypes(1).sessionIndex == 1, :);
mid_dFF_cued = byTypes(1).phRaster(byTypes(1).sessionIndex == 2, :);
fourth_dFF_cued = byTypes(1).phRaster(byTypes(1).sessionIndex == 4, :);
late_dFF_cued = byTypes(1).phRaster(byTypes(1).sessionIndex == 5, :);
late_dFF_uncued = byTypes(5).phRaster(byTypes(5).sessionIndex == 5, :);
ensureFigure('mid_singleTrials', 1); plot(xdata, mid_dFF_cued(52:61, :)', 'color', [0.6 0.6 0.6]); % not gonna get much better
hold on;
plot(xdata, mean(mid_dFF_cued), 'k');
formatFigure;
saveas(gcf, 'singleTrials.fig');
saveas(gcf, 'singleTrials.epsc');

%%
ensureFigure('cued_ex_phRaster', 1);
imagesc(xdata, 1:size(fourth_dFF_cued, 1), fourth_dFF_cued, [-.025 .025]);
set(gca, 'TickDir', 'out');
formatFigure;
saveas(gcf, 'phRaster.fig');
saveas(gcf, 'phRaster.epsc');
%%
ensureFigure('early_cued_phRaster', 1);
imagesc(xdata, 1:size(early_dFF_cued, 1), early_dFF_cued, [-.025 .025]);
set(gca, 'TickDir', 'out');
formatFigure;
saveas(gcf, 'early_phRaster.fig');
saveas(gcf, 'early_phRaster.epsc');
%%
ensureFigure('expectency', 1);
mean_cued = mean(late_dFF_cued);
sem_cued = std(late_dFF_cued) / sqrt(size(late_dFF_cued, 1));
mean_uncued = mean(late_dFF_uncued);
sem_uncued = std(late_dFF_uncued) / sqrt(size(late_dFF_uncued, 1));
boundedline(xdata, mean_cued, sem_cued, 'b', xdata, mean_uncued, sem_uncued, 'c');
formatFigure;
saveas(gcf, 'expectency.fig');
saveas(gcf, 'expectency.epsc');
%%
firstPoint = 10;
xDataNew = xdata(firstPoint:end);
ylim = [-0.005 0.020]; xlim = [-3 6];
figSize = [4 2];
ah = zeros(3,1);
% early 
saveName = 'early_phHist';
ensureFigure(saveName, 1); ah(1) = axes;
mean_cued = mean(early_dFF_cued);
sem_cued = std(early_dFF_cued) / sqrt(size(early_dFF_cued, 1));
[lines, patches] = boundedline(xDataNew, mean_cued(firstPoint:end), sem_cued(firstPoint:end), 'cmap', [171 55 214]/256);
set(lines, 'LineWidth', 2);
set(gca, 'YLim', ylim, 'XLim', xlim, 'YTickLabel', '');
formatFigurePoster(figSize);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    
% middle
saveName = 'mid_phHist';
ensureFigure(saveName, 1); ah(2) = axes;
mean_cued = mean(mid_dFF_cued);
sem_cued = std(mid_dFF_cued) / sqrt(size(mid_dFF_cued, 1));
[lines, patches] = boundedline(xDataNew, mean_cued(firstPoint:end), sem_cued(firstPoint:end), 'cmap', [171 55 214]/256);
set(lines, 'LineWidth', 2);
set(gca, 'YLim', ylim, 'XLim', xlim, 'YTickLabel', '');
formatFigurePoster(figSize);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    
% late
saveName = 'late_phHist';
ensureFigure(saveName, 1);
mean_cued = mean(late_dFF_cued); ah(3) = axes;
sem_cued = std(late_dFF_cued) / sqrt(size(late_dFF_cued, 1));
[lines, patches] = boundedline(xDataNew, mean_cued(firstPoint:end), sem_cued(firstPoint:end), 'cmap', [171 55 214]/256);
set(lines, 'LineWidth', 2);
set(gca, 'YLim', ylim, 'XLim', xlim, 'YTickLabel', '');
formatFigurePoster(figSize);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    





%% SFN 2017
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SFN_2017\PE_Hypothesis';
clim = [-0.02 0.015];
ensureFigure('early_cued_phRaster', 1);
colormap('parula');
imagesc(xdata, 1:size(early_dFF_cued, 1), early_dFF_cued, [clim]);
set(gca, 'TickDir', 'out', 'Box', 'off', 'XTick', [], 'YTick', []);
formatFigurePoster([4,3]);
saveas(gcf, fullfile(savepath, 'early_phRaster'), 'fig');
saveas(gcf, fullfile(savepath, 'early_phRaster'), 'jpeg');
saveas(gcf, fullfile(savepath, 'early_phRaster'), 'epsc');


ensureFigure('mid_cued_phRaster', 1);
imagesc(xdata, 1:size(mid_dFF_cued, 1), mid_dFF_cued, [clim]);
set(gca, 'TickDir', 'out', 'Box', 'off', 'XTick', [], 'YTick', []);
set(gca, 'YLim', [10 80]);
formatFigurePoster([4,3]);
saveas(gcf, fullfile(savepath, 'mid_phRaster'), 'fig');
saveas(gcf, fullfile(savepath, 'mid_phRaster'), 'jpeg');
saveas(gcf, fullfile(savepath, 'mid_phRaster'), 'epsc');

ensureFigure('late_cued_phRaster', 1);
imagesc(xdata, 1:size(late_dFF_cued, 1), late_dFF_cued, [clim]);
set(gca, 'TickDir', 'out', 'Box', 'off', 'XTick', [], 'YTick', []);
formatFigurePoster([4,3]);
saveas(gcf, fullfile(savepath, 'late_phRaster'), 'fig');
saveas(gcf, fullfile(savepath, 'late_phRaster'), 'jpeg');
saveas(gcf, fullfile(savepath, 'late_phRaster'), 'epsc');

%%
save('McKnight.mat', 'allTE', 'TE', 'byTypes');
% save('McKnight_sessions.mat', 'sessions', '-v7.3');

%%
    return
    cd('C:\Users\Adam\Dropbox\KepecsLab\_Fitz\McKnight_Poster');
    ensureFigure('lickHist_McKnight', 1);
    axes;
    bpLickHist(session.SessionData, [1 5], [1 1], [-2 6 0.25],...
        'DeliverStimulus', 'PreCsRecording', 'PostUsRecording', {'b', 'c'}, [], gca);
    
    xlabel('time (s)'); ylabel('Lick rate (1/s)');
    formatFigure;
    saveas(gcf, 'lickHist_McKnight.fig');
    saveas(gcf, 'lickHist_McKnight.epsc');    