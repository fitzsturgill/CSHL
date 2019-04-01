






% %% SIGNAL CORRELATIONS: show that cue and reward responses are correlated on a trial-by-trial basis between let and right BLA
% 
% saveName = 'LeftRight_BLA_SIGNAL_correlations';
% ensureFigure(saveName, 1);
% hitTrials = TE.licks_cs.rate > 0;
% linecolors = [0 0 1; 0 1 1; 1 0 0; 0 1 0; mycolors('shock')];         
% trialSets = [Odor2Valve1Trials & hitTrials & rewardTrials, uncuedReward, Odor2Valve2Trials & punishTrials, Odor2Valve2Trials & shockTrials];
% allTrials = sum(trialSets, 2) ~= 0;
% % xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];
% 
% subplot(1,2,1); hold on; subplot(1,2,2); hold on;
% for counter = 1:size(trialSets, 2)
%     h=[];
%     h(1) = subplot(1,2,1); %set(gca, 'XLim', xlims(1,:));
% %     allTrials = allTrials | trialSets{counter};
%     xData = TE.phPeakMean_cs(2).data(trialSets(:, counter)); yData = TE.phPeakMean_cs(1).data(trialSets(:, counter)); 
% %     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
%     scatter(xData, yData, 8, linecolors(counter, :), '.');
%     % fit for cs
%     fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
%     fob = fit(xData, yData, 'poly1', fo); 
%     fph=plot(fob); legend off; %,'predfunc'); legend off;
%     set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
% 
%     h(2) = subplot(1,2,2); %set(gca, 'XLim', xlims(2,:));
%     xData = TE.phPeakMean_us(2).data(trialSets(:, counter)); yData = TE.phPeakMean_us(1).data(trialSets(:, counter)); 
% %     scatter(TE.phPeakMean_us(2).data(trialSets{counter}), TE.phPeakMean_us(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');    
%     scatter(xData, yData, 8, linecolors(counter, :), '.');
%     % fit for us
%     fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
%     fob = fit(xData, yData, 'poly1', fo); 
%     fph=plot(fob); legend off;% ,'predfunc'); legend off;
%     set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
%     sameXYScale(h);
% end
% subplot(1,2,1);
% textBox('Cue', [], [0.5 0.95], 8);
% xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
% ylabel('\fontsize{8}Right (\fontsize{12}\sigma\fontsize{8}-baseline)');
% subplot(1,2,2);
% textBox('Outcome', [], [0.5 0.95], 8);
% xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
% ylabel('');
% 
% formatFigurePublish('size', [2.5 1.1]);
% 
% if saveOn 
%     export_fig(fullfile(savepath, saveName), '-eps');
% end
% 
% %% NOISE CORRELATIONS: show that cue and reward responses are correlated on a trial-by-trial basis between let and right BLA
% 
% saveName = 'LeftRight_BLA_NOISE_correlations';
% ensureFigure(saveName, 1);
% hitTrials = TE.licks_cs.rate > 0;
% linecolors = [0 0 1; 0 1 1; 1 0 0; 0 1 0; mycolors('shock')];         
% trialSets = [Odor2Valve1Trials & hitTrials & rewardTrials, uncuedReward, Odor2Valve2Trials & punishTrials, Odor2Valve2Trials & shockTrials];
% allTrials = sum(trialSets, 2) ~= 0;
% % xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];
% 
% subplot(1,2,1); hold on; subplot(1,2,2); hold on;
% for counter = 1:size(trialSets, 2)
%     h=[];
%     h(1) = subplot(1,2,1); %set(gca, 'XLim', xlims(1,:));
% %     allTrials = allTrials | trialSets{counter};
%     xData = TE.phPeakMean_cs(2).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_cs(2).data(trialSets(:, counter))); 
%     yData = TE.phPeakMean_cs(1).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_cs(1).data(trialSets(:, counter))); 
% %     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
%     scatter(xData, yData, 8, linecolors(counter, :), '.');
%     % fit for cs
%     fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
%     fob = fit(xData, yData, 'poly1', fo); 
%     fph=plot(fob); legend off; %,'predfunc'); legend off;
%     set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
% 
%     h(2) = subplot(1,2,2); %set(gca, 'XLim', xlims(2,:));
%     xData = TE.phPeakMean_us(2).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_us(2).data(trialSets(:, counter))); 
%     yData = TE.phPeakMean_us(1).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_us(1).data(trialSets(:, counter))); 
% %     scatter(TE.phPeakMean_us(2).data(trialSets{counter}), TE.phPeakMean_us(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');    
%     scatter(xData, yData, 8, linecolors(counter, :), '.');
%     % fit for us
%     fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
%     fob = fit(xData, yData, 'poly1', fo); 
%     fph=plot(fob); legend off;% ,'predfunc'); legend off;
%     set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
%     sameXYScale(h);
% end
% subplot(1,2,1);
% textBox('Cue', [], [0.5 0.95], 8);
% xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
% ylabel('\fontsize{8}Right (\fontsize{12}\sigma\fontsize{8}-baseline)');
% subplot(1,2,2);
% textBox('Outcome', [], [0.5 0.95], 8);
% xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
% ylabel('');
% 
% formatFigurePublish('size', [2.5 1.1]);
% 
% if saveOn 
%     export_fig(fullfile(savepath, saveName), '-eps');
% end