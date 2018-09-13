% compile photometry averages, make grand average
fdField = 'ZS';

ensureFigure('test', 1);
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadExperiment(DB, animal);
    if ~success
        disp('wtf');
        continue
    end    
    disp(animal);
    for channel = 1:2
        avgData = phAverageFromTE(TE, {rewardTrials & csPlusTrials & hitTrials, uncuedReward}, channel,...
            'FluorDataField', fdField, 'window', [1, 7]); %high value, reward
            grand(channel).reward.components(:,:,counter) = avgData.Avg;
            grand(channel).reward.xData = avgData.xData(1,:); 
        if channel == 1
            subplot(3,3,counter)
            plot(avgData.xData(1,:), avgData.Avg');
        end




        avgData = phAverageFromTE(TE, {csPlusTrials & hitTrials, csPlusTrials & missTrials, csMinusTrials & CRTrials}, channel,...
            'FluorDataField', fdField, 'window', [-3, 7]); %high value, reward
            grand(channel).cue.components(:,:,counter) = avgData.Avg;
            grand(channel).cue.xData = avgData.xData(1,:);         
    end
end