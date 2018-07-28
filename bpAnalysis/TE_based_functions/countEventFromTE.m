function eventCount = countEventFromTE(TE, event, window, zeroTimes, varargin)
    % used in conjunction with bpAddEventAsTrialEvent
    % window- can be either a 1x2 row vector or a nTrials x 2 matrix (to
    % provide specific windows for every trial)


    defaults = {...    
        'referenceFromEnd', 0;... % if you provide zeroTimes as a cell array of times (say extracted from state times), you can use start of first or end of last
        };    
    [s, ~] = parse_args(defaults, varargin{:});
    %%

    nTrials = length(TE.trialNumber);
    
    eventCount = struct(...
        'count', zeros(nTrials,1),...
        'rate', zeros(nTrials,1),...
        'duration', zeros(nTrials,1),...
        'settings', s...
        );
    %%
    if isscalar(zeroTimes) % for example if you use bpod time directly you don't need a zero (zero is start of bpod time)
        zeroTimes = repmat(zeroTimes, nTrials, 1);
    end
    
    if iscell(zeroTimes)
        if ~s.referenceFromEnd
            zeroTimes2 = cellfun(@(x) x(1), zeroTimes); 
        else
            zeroTimes2 = cellfun(@(x) x(end), zeroTimes);
        end
    else
        zeroTimes2 = zeroTimes; % not yet tested
    end
    
    if size(window, 1) == 1
        window = repmat(window, nTrials, 1);
    end
    
    %%
    for trial = 1:nTrials
        trialEvents = TE.(event){trial};
        if ~isnan(trialEvents(1))
            trialWindow = window(trial,:) + zeroTimes2(trial);
            eventCount.count(trial) = length(find(trialWindow(1) < trialEvents & trialEvents < trialWindow(2)));
            eventCount.duration(trial) = diff(trialWindow);
            eventCount.rate(trial) = eventCount.count(trial) / eventCount.duration(trial);
        end
    end
            
    
    