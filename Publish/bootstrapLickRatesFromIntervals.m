function resampled = bootstrapLickRatesFromIntervals(intervals, durations, nBoot)
%{ 
estimate sample distribution of lick rates across windows
of specified duration by resmpling inter-lick intervals
Inputs: 
intervals- nReverals x nTrials cell array of interlick intervals (sum <=
duration), if empty means there are no licks at all
durations- durations of windows for a given trial
Outputs:
resampled = nReversals x nTrials x nBootstraps
%}
      

nRevs = size(intervals, 1);
nTrials = size(intervals, 2);
resampled = NaN(nRevs, nTrials, nBoot);
h = waitbar(0, 'bootstrapping lick rates');
for revCounter = 1:nRevs
    for trialCounter = 1:nTrials
        theseIntervals = intervals{revCounter, trialCounter};
        if isnan(theseIntervals)
            continue
        elseif isempty(theseIntervals) % no intervals means no events
            resampled(revCounter, trialCounter, :) = 0;
            continue
        end
        thisDuration = durations(revCounter, trialCounter);
        nIntervals = length(theseIntervals);
        for bootCounter = 1:nBoot
            eventSum = 0;
            nEvents = 0;
            while 1
                eventSum = eventSum + theseIntervals(randi(nIntervals, 1));
                if eventSum >=  thisDuration
                    break
                end
                nEvents = nEvents + 1;
            end
            resampled(revCounter, trialCounter, bootCounter) = nEvents / thisDuration;
        end
    end
    waitbar(revCounter / nRevs);
end
close(h);
        