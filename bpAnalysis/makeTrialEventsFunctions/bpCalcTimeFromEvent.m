function timeMatrix = bpCalcTimeFromEvent(TE, event, varargin)

% process wheel running speed measurements from wheel connected to rotary
% encoder.  Assumes that rotary pulses are unidirectional.
    %% optional parameters, first set defaults
    defaults = {...
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (store data in cell array)
        'trialStartField', 'trialStartTimeStamp';...
        'dataStartField', 'Start';...
        'Fs', 20;...
        'duration', 30;...
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    nTrials = length(TE.filename);
    nSamples = s.duration * s.Fs;
    %% generate absolute time vectors for each trial
    timeMatrix = repmat(linspace(0, s.duration - 1/s.Fs, nSamples), nTrials, 1);
    timeMatrix = bsxfun(@plus, timeMatrix, TE.(s.trialStartField));
    timeMatrix = bsxfun(@plus, timeMatrix, TE.(s.dataStartField));
    %% calculate absolute event times
    eventAbs = TE.(s.event);    
    for counter = 1:nTrials
        trailStart = TE.(s.trialStartField)(counter);
        eventAbs{counter} = eventAbs{counter} + trialStart;        
    end
    
    %% generate event time annotated matrix to be subtracted from time matrix
    annotatedMatrix = zeros(nTrials, nSamples);    
    for counter = 1:nTrials
        eventTimes = eventAbs{counter}(:,1);
        timeVector = timeMatrix(counter, :);
        annotatedVector = zeros(size(timeVector));
        eventTimesWithinTrial = eventTimes(eventTimes >= min(timeVector) & eventTimes < max(timeVector)); % only events occuring within time span of data for a given trial
        if counter == 1
            eventLastTrial = inf;
        else
            eventLastTrial = eventAbs{counter - 1}(end,1);
        end
        nEvents = size(eventTimesWithinTrial, 1);
        eventIndices = zeros(nEvents, 1);
        for i = 1:nEvents
            eventIndices(i) = nearest(timeVector, eventTimesWithinTrial(i));
        end
        % determine last event preceding data in a given trial
        precedingEvents = eventTimes(eventTimes < min(timeVector));
        if isempty(precedingEvents)
            precedingEvent = eventLastTrial;
        else
            precedingEvent = precedingEvents(end);
        end
        
        if eventIndices(1) > 1
            annotatedVector(1:eventIndices(1) - 1) = precedingEvent;
        end
        
        if length(eventIndices > 1)
            for j = 1:length(eventIndices) - 1
                annotatedVector(eventIndices(j):eventIndices(j+1) - 1) = eventTimesWithinTrial(j);
            end
        end
        
        annotatedVector(eventIndices(end):end) = eventTimesWithinTrial(end);
        
        annotatedMatrix(counter, :) = annotatedVector;
    end
        
        
        
    