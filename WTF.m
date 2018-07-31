% ensureFigure('test', 1);
% plot(TE.Photometry.bleachFit(counter, channel).templateX, TE.Photometry.bleachFit(counter, channel).trialTemplate, 'g'); hold on;      
% plot(TE.Photometry.bleachFit(counter, channel).fitX, TE.Photometry.bleachFit(counter, channel).trialTemplateFull, 'b'); 

% [xData avgData] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, neutralTrials, uncuedReward, uncuedPunish}, 1,...
%             'FluorDataField', 'ZS', 'window', [1, 7]);
        
        
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, 'method', 'mean', 'phField', 'ZS');