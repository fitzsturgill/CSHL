function counts = histCountsByTrial(times, trials, binEdges)
% '2d' histcounts where nRows = number of unique trials
% assumes consistent trial timing (~uniformOutput = 1)
    validTrials = unique(trials);
    nValidTrials = length(validTrials);
    
    counts = zeros(nValidTrials, length(binEdges) - 1);
    
    for counter = 1:nValidTrials
        timesThisTrial = times(trials == validTrials(counter));
%         if any(isnan(timesThisTrial))
%             disp('weird');
%         end
        counts(counter, :) = histcounts(timesThisTrial, binEdges);
    end