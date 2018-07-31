% ensureFigure('test', 1);
% plot(TE.Photometry.bleachFit(counter, channel).templateX, TE.Photometry.bleachFit(counter, channel).trialTemplate, 'g'); hold on;      
% plot(TE.Photometry.bleachFit(counter, channel).fitX, TE.Photometry.bleachFit(counter, channel).trialTemplateFull, 'b'); 

% [xData avgData] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, neutralTrials, uncuedReward, uncuedPunish}, 1,...
%             'FluorDataField', 'ZS', 'window', [1, 7]);
        
        
phRasterFromTE(TE, csPlusTrials & (day1Trials | day2Trials) & rewardTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');