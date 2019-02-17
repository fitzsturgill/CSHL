function [eventTimes, eventTrials] = extractEventTimesFromTE(TE, trials, event, varargin)
% eventTimes, size nEvents x 1 vector of event times relative to zero as defined by
% zeroField or zeroTimes arguments
% eventTrials, size nEvents x 1 vector of trial numbers for events
%     8/2016
%     7/2018, updated to use zeroTimes and window (optionally instead of
%     start, end, and zero fields)

    defaults = {...
        'startField', [];...
        'zeroField', [];... % required
        'endField', [];...
        'zeroTimes', [];... % by convention, of size nTotalTrials (from TE)  x 1   (not just length of trials argument)
        'window', [];... % either of size 1 x 2 or nTotalTrials x 2   
        'trialNumbering', 'global';... % 'singleSession', 'consecutive' (doesn't repeat numbers when switching sessions), 'global' cross sessions and cross-conditions
        };
    
    [s, ~] = parse_args(defaults, varargin{:});
    assert(~isempty(s.zeroField) || ~isempty(s.zeroTimes), 'must supply either zeroField or zeroTimes');
    
    eventTimes = [];
    eventTrials = [];

    if islogical(trials)
        trials = find(trials);
    end
    
    if size(s.window, 1) == 1
        s.window = repmat(s.window, length(TE.filename), 1);
    end
    
    for counter = 1:length(trials)
        trial = trials(counter);
        
        if ~isempty(s.zeroField)
            zeroTime = TE.(s.zeroField){trial}(1);        
        elseif iscell(s.zeroTimes)
            zeroTime = s.zeroTimes{trial}(1);
        else
            zeroTime = s.zeroTimes(trial);
        end
        
        if ~isempty(s.startField)
            startTime = TE.(s.startField){trial}(1);
        elseif ~isempty(s.window)
            startTime = zeroTime + s.window(trial, 1);
        else
            startTime = -Inf;
        end
        
        if ~isempty(s.endField)
            endTime = TE.(s.endField){trial}(end);
        elseif ~isempty(s.window)
            endTime = zeroTime + s.window(trial, 2);            
        else
            endTime = Inf;
        end
                      
        trialEventTimes = TE.(event){trial};
        trialEventTimes = trialEventTimes((startTime <= trialEventTimes) & (trialEventTimes < endTime));
        if isempty(trialEventTimes)
            trialEventTimes = NaN;
        end
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
    

        