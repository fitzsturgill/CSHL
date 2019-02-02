saveName = 'Aversive_blinkAverages';
ensureFigure(saveName, 1);
xData = TE.eyeAvg.xData;
blinkData = TE.eyeAvg.EyeAvgNorm;
blinkData(isnan(blinkData)) = 0;
axes; hold on; grid on;
boundedline(xData(:), [mean(blinkData(trialsByType{3}, :)); mean(blinkData(trialsByType{4}, :)); mean(blinkData(trialsByType{6}, :))]',...
    permute([std(blinkData(trialsByType{3}, :)) ./ sqrt(sum(trialsByType{3})); std(blinkData(trialsByType{4}, :)) ./ sqrt(sum(trialsByType{3})); std(blinkData(trialsByType{6}, :)) ./ sqrt(sum(trialsByType{3}))], [2 3 1]));