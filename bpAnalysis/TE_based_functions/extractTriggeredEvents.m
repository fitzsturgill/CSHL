function tre = extractTriggeredEvents(TE, event, TS, varargin)

% extracts times relative to a time stamps and trialNumbers of an event derived from a
% TE structure 

% TS, time stamps, cell array nTrials x 1, assumes that time stamps are in
% the first column of cell array contents for a given trial (so that you
% can use output of bpAddStateAsTrialEvent directly)
% window, time window around time stamp to extract

    defaults = {...
        'trials', [];... % if empty, sift through all trials
        'trialNumbering', 'global';... % 'singleSession', 'consecutive' (doesn't repeat numbers when switching sessions), 'global' cross sessions and cross-conditions
        'window', [-2 2];...
        };
    
    [s, ~] = parse_args(defaults, varargin{:});
    
    if isempty(s.trials)
        s.trials = 1:length(TE.filename);
    end

    % find trials that contain a time stamp
    validTrials = cellfun(@(x) ~isnan(x(:,1)), TS, 'UniformOutput', 0); % use time stamps in first column    
    validTrials = cellfun(@(x) sum(x), validTrials);
%     nStamps = sum(validTrials(s.trials));
    if islogical(s.trials)
        s.trials = find(s.trials);
    end
    validTrials = intersect(find(validTrials), s.trials);
    validTrials = validTrials(:)'; % ensure row vector
    tre.eventTimes = [];
    tre.eventTrials = [];
    tre.eventStamps = [];
    tre.eventStampIndices = [];    
    stampCounter = 1;
    for trial = validTrials
        stamps = (TS{trial}(:,1));  % use first column of time stamps (in case you are using output of bpAddStateAsTrialEvent);
        stamps = stamps(:)';
        for stamp = stamps
            trialWindow = s.window + stamp;
            trialEventTimes = TE.(event){trial};
            trialEventTimes = trialEventTimes((trialWindow(1) <= trialEventTimes) & (trialEventTimes < trialWindow(2)));
            trialEventTimesZeroed = trialEventTimes - stamp;
            nEvents = numel(trialEventTimesZeroed);

            trialEventTimesZeroed = trialEventTimesZeroed(:); % ensure column
            tre.eventTimes = [tre.eventTimes; trialEventTimesZeroed];
            tre.eventTrials = [tre.eventTrials; repmat(trial, nEvents, 1)];
            tre.eventStampIndices = [tre.eventStampIndices; repmat(stampCounter, nEvents, 1)];
            tre.eventStamps = [tre.eventStamps; repmat(stamp, nEvents, 1)];            
            stampCounter = stampCounter + 1;
        end
    end
            

    
    
    