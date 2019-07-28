%% Development code below: normalize responses to punishment by those to reward and compare between left and right BLA


FluorField_us = 'phPeakMean_us';

saveName = sprintf('LeftRight_PunNorm_%s', FluorField_us);
ensureFigure(saveName, 1);
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
    linecolors = [0 0 1; 1 0 0; mycolors('shock')];         
    trialSets = {'rew', 'puff', 'shock'};
    trialSetNames = {'Reward', 'Air Puff', 'Shock'};
    allTrials = sum(trialSets, 2) ~= 0;    
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

        h(end + 1) = scatter(xData, yData, 20, linecolors(counter, :), '.', 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);

        xDataStruct(counter).Mean  = nanmean(xData);
        xDataStruct(counter).SEM = nanstd(xData) / sqrt(sum(isfinite(xData)));
        yDataStruct(counter).Mean = nanmean(yData);
        yDataStruct(counter).SEM = nanstd(yData) / sqrt(sum(isfinite(yData)));
    end
    addUnityLine(gca, [0.3 0.3 0.3]);

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
        errorbar(xDataStruct(counter).Mean, yDataStruct(counter).Mean, -yDataStruct(counter).SEM, yDataStruct(counter).SEM,-xDataStruct(counter).SEM, xDataStruct(counter).SEM, 'Color', linecolors(counter, :), 'LineWidth', 1, 'CapSize', 2, 'Marker', 'none');
        h(end + 1) = plot(xDataStruct(counter).Mean, yDataStruct(counter).Mean, '-', 'Color', linecolors(counter, :));
    end
    
    if acounter == 2
        set(gca, 'YLim', [-5 6], 'XLim', [-5 6]);
    end
    sameXYScale(gca);
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