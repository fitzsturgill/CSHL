%% omitReward, aligned by next lick, if present
window = [-5 2];
nSamples = diff(window) * 20;
nSites = length(DB.animals) * 2;

omitReward_lickAligned = struct(...
    'avgData', NaN(nSites, nSamples),...
    'semData', NaN(nSites, nSamples),...
    'xData', [],...
    'animal', [],...
    'ch', []...
    );


for counter = 1:length(DB.animals)
    animal = DB.animals{counter};

    dbLoadAnimal(DB, animal);
    
    % omitReward
    trials = neutralTrials & Odor2Valve1Trials & ismember(TE.BlockNumber, [2 3]);          
    for ch = 1:length(TE.Photometry.data)
        thisSite = (counter - 1)*2 + ch;
        avgData = phAverageFromTE(TE, trials, ch, 'zeroTimes', cellfun(@(x) x(1), TE.Outcome) + TE.lickLatency_us, 'window', window, 'FluorDataField', 'ZS');
        omitReward_lickAligned.avgData(thisSite,:) = avgData.Avg;
        omitReward_lickAligned.semData(thisSite,:) = avgData.SEM;
        if isempty(omitReward_lickAligned.xData)
            omitReward_lickAligned.xData = avgData.xData;
        end
    end
end