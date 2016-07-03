function [lickRates, lickRatesN, binCenters] =  bpLickCounts(SessionData, type, outcome, zeroField, binSpecs, startField, endField)
    
    %Arguments:
    % binSpecs 3 element Vector of form [start, end, width]
    
   % example [lickRates, lickRatesN, binCenters] =  bpLickCounts(SessionData, 1, 2, 'StartStimulus', [-4, 4, 0.5], 'PreTrialRecording', 'PostTrialRecording')
   % !!!!! zeroField, startField and endField-   can also take the form of
   % cell arrays in order to indicate whether to use the FIRST or LAST
   % event, e.g. startField = {'restartNoLick','first'}; the keywords are
   % 'first' and 'last' for this (second element in cell array)


    trials = bpFilterTrials(SessionData, type, outcome);
    binEdges = linspace(binSpecs(1), binSpecs(2), round((binSpecs(2) - binSpecs(1))/binSpecs(3))); % if interval isn't divisible by width then adjust number of bins
    binSpecs(3) = binEdges(2) - binEdges(1); % and adjust width if necessary
    % initialize 
    lickRates = zeros(length(trials), length(binEdges) - 1);
    lickRatesN = zeros(1, length(binEdges) - 1);
    binCenters = binEdges(1:end-1) + binSpecs(3)/2;
    for i = 1:length(trials)
        trial = trials(i);
        if ~iscell(zeroField)
            zeroTime = SessionData.RawEvents.Trial{trial}.States.(zeroField)(1,1); % always use first instance and start of zeroField for zeroing (each state has a start and end timestamp)
        else
            zeroTime = getTime(zeroField); % if zeroField provided as a cell, you are using the 'first' or 'last' keywords to indicate which instance of that event to use for zeroing
        end

        % determine which bins are actually included in this trial so that
        % you don't depress lick rates without justification due to
        % relatively short trials
        if ~iscell(startField)
            start = SessionData.RawEvents.Trial{trial}.States.(startField)(1,1) - zeroTime;
        else
            start = getTime(startField) - zeroTime;
        end
        if ~iscell(endField)
            finish = SessionData.RawEvents.Trial{trial}.States.(endField)(end,2) - zeroTime;
        else
            finish = getTime(endField) - zeroTime;
        end            
        whichBins = binEdges > start & binEdges <= finish;
        trialBins = binEdges(whichBins);        
        lickRatesN = lickRatesN + whichBins(1:end-1);        
        if isfield(SessionData.RawEvents.Trial{trial}.Events, 'Port1In')        
            theseLicks = SessionData.RawEvents.Trial{trial}.Events.Port1In;
            theseLicks = theseLicks - zeroTime;
            [trialRates, ~] = histcounts(theseLicks, trialBins);
            trialRates = trialRates ./ binSpecs(3);
            trialRates = [zeros(1, find(whichBins, 1) - 1) trialRates];
            if length(trialRates) < length(lickRatesN)
                trialRates = [trialRates zeros(1, length(lickRatesN) - length(trialRates))]; % make up difference if a short trial
            end
            lickRates(i, :) = trialRates;
        else
            lickRates(i, :) = 0;
        end
    end

    function out = getTime(timeCell)
        whichEvent = timeCell{2};
        eventType = timeCell{1};
        switch whichEvent
            case 'first'
                out = SessionData.RawEvents.Trial{trial}.States.(eventType)(1,1); % first event, beginning of
            case 'last'
                out = SessionData.RawEvents.Trial{trial}.States.(eventType)(end,2); % last event, end of
        end
    end
end