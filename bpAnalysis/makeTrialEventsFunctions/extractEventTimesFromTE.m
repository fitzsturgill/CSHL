function [eventTimes, eventTrials] = extractEventTimesFromTE(TE, trials, event, varargin)
    % first draft 8/2016
    defaults = {...
        'startField', [];...
        'zeroField', [];... % required
        'endField', [];...
        'trialNumbering', 'global';... % 'singleSession', 'consecutive' (doesn't repeat numbers when switching sessions), 'global' cross sessions and cross-conditions
        };
    
    [s, ~] = parse_args(defaults, varargin{:});
    if isempty(s.zeroField)
        error('zeroField must be supplied');
    end
    eventTimes = [];
    eventTrials = [];

    if islogical(trials)
        trials = find(trials);
    end
    for counter = 1:length(trials)
        trial = trials(counter);
        startTime = TE.(s.startField){trial}(1);
        zeroTime = TE.(s.zeroField){trial}(1);
        endTime = TE.(s.endField){trial}(end);        
        trialEventTimes = TE.(event){trial};
        trialEventTimes = trialEventTimes((startTime < trialEventTimes) & (trialEventTimes < endTime));
        trialEventTimesZeroed = trialEventTimes - zeroTime; 
        if size(trialEventTimesZeroed, 2) > 1
            trialEventTimesZeroed = trialEventTimesZeroed';
        end
        eventTimes = [eventTimes; trialEventTimesZeroed];
        switch s.trialNumbering
            case 'global'
                eventTrials = [eventTrials; repmat(trial, length(trialEventTimesZeroed), 1)];
            case 'consecutive'
                eventTrials = [eventTrials; repmat(counter, length(trialEventTimesZeroed), 1)];
            case 'singleSession'
                eventTrials = [eventTrials; repmat(TE.trialNumber(trial), length(trialEventTimesZeroed), 1)];
        end
    end
        