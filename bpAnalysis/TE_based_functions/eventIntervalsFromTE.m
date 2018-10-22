function intervals = eventIntervalsFromTE(TE, event, window, zeroTimes, varargin)
    % used in conjunction with bpAddEventAsTrialEvent
    % window- can be either a 1x2 row vector or a nTrials x 2 matrix (to
    % provide specific windows for every trial)

    % Important- you can omit zeroTimes if window is specified in Bpod time
    % (i.e. relative to Bpod trial start)
    defaults = {...    
        'referenceFromEnd', 0;... % if you provide zeroTimes as a cell array of times (say extracted from state times), you can use start of first or end of last
        };    
    [s, ~] = parse_args(defaults, varargin{:});
    %%

    nTrials = length(TE.trialNumber);
    
    intervals = struct(...
        'rate', zeros(nTrials,1),...
        'count', zeros(nTrials, 1),...
        'duration', zeros(nTrials,1),...
        'settings', s...
        );
    intervals.intervals = cell(nTrials, 1);
    %%
    if nargin < 4 % if you use bpod time directly you don't need a zero (zero is start of bpod time)
        zeroTimes = zeros(nTrials, 1);
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
            eventIndices = find(trialWindow(1) < trialEvents & trialEvents < trialWindow(2));
            intervals.count(trial) = length(eventIndices);
            intervals.intervals{trial} = diff([trialWindow(1) trialEvents(eventIndices) trialWindow(2)]);
            intervals.duration(trial) = diff(trialWindow);
            intervals.rate(trial) = intervals.count(trial) / intervals.duration(trial);
        end
    end
            
    
    