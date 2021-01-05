saveOn = 1;
savePath = 'C:\Users\Fitz\Dropbox\KepecsLab\_Fitz\Manuscripts\ChAT_paper\Figures\Cell\data_balazs_panna\';
figSize = [8 8];
%% cue
open(fullfile(savePath, 'CuedOutcome_ScatterPlot_Cue.fig'));
saveName = 'CuedOutcome_ScatterPlot_Cue_final';
ebData = get(gca, 'Children');
% minmax = [min(min(ebData.YData)) max(max(ebData.YData))];
axlims = [0.5 25];
set(gca, 'YLim', axlims, 'XLim', axlims);
set(gca, 'XScale', 'log', 'YScale', 'log');
addUnityLine;
formatFigurePublish('size', figSize);
set(gca, 'YTick', [1 10], 'XTick', [1 10]);
if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));   
end  

%% scrape data and make whisker plot to match existing ones for reward and punish
cueDiff = ebData.YData - ebData.XData;
saveName = 'cue_boxPlot';
ensureFigure(saveName);

boxplot(cueDiff, 'BoxStyle', 'filled', 'Colors', [0.5 0.5 0.5], 'Symbol', '');
set(gca, 'YLim', [-1 3], 'YTick', [0 2]);
formatFigurePublish('size', [1 2]);
if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));   
end  
%% behavior example
open(fullfile(savePath, 'behavior_example.fig'));
saveName = 'behavior_example_final';
axs = get(gcf, 'Children');
axs(2:4).delete;
axs = axs([1 5 6]);
%% 
axes(axs(1));
set(gca, 'XTick', [-3 0 3]);
xlabel(''); ylabel('');
set(axs(2:3), 'YAxisLocation', 'left');
formatFigurePublish('size', [2 3]);
if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));   
end  


%% Outcome Valence
open(fullfile(savePath, 'CuedOutcome_ScatterPlot_Reward.fig'));
saveNameR = 'CuedOutcome_ScatterPlot_Reward_final';
set(gcf, 'Name', saveNameR);
fhr = gcf;

open(fullfile(savePath, 'CuedOutcome_ScatterPlot_Punishment.fig'));
saveNameP = 'CuedOutcome_ScatterPlot_Punishment_final';
set(gcf, 'Name', saveNameP);
fhp = gcf;

% get outcome data for reward and punishment (from most common condition)
figure(fhr);
figdata = get(gca, 'Children');
reward_mean = figdata.XData;
reward_sem = figdata.XPositiveDelta;

figure(fhp);
figdata = get(gca, 'Children');
punish_mean = figdata.YData;
punish_sem = figdata.YPositiveDelta;

saveName = 'CuedOutcome_outcomeValence_spikes';
ensureFigure(saveName); 
axes; hold on;
scatter(reward_mean,punish_mean, 36, [0 0 0], 'filled');
errorbar(reward_mean,punish_mean, -punish_sem, punish_sem,...
    -reward_sem, reward_sem, '.', 'Color', [0 0 0]);
scatter(reward_mean,punish_mean, 36, [0.5 0.5 0.5], 'filled');
set(gca, 'YLim', [-10 40], 'XLim', [-10 40], 'YTick', [0 20 40], 'XTick', [0 20 40]);
addOrginLines;
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));   
end  

