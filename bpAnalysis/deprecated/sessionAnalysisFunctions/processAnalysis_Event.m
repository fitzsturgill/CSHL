function analysis = processAnalysis_Event(SessionData, analysis, event, varargin)
    % collates trial events by trial type and trial outcome
    % if you want to zero events multiple times, 
    % output them into different analysis structures 
    % ex. session.analysis_zero_by_stim = struct();
    % session.analysis_zero_by_stim =
    % processAnalysis_Event(session.SessionData, session.analysis_zero_by_stim, 'Port1In', 'zeroField',
    % 'DeliverStimulus');
    
    % 
    if ~isstruct(analysis)
        disp('error in processAnalysis_Event');
        analysis = [];
    end
    
    analysis.(event) = struct(...
        'data', [],... % 
        'settings', {varargin}...
        );
    
    %default values
    zeroField = '';
    startField = '';
    endField = '';
    trialTypes = unique(SessionData.TrialTypes);
    outcomes = unique(SessionData.TrialOutcome);    

    % parse input parameter pairs
    counter = 1;
    while counter+1 <= length(varargin) 
        prop = varargin{counter};
        val = varargin{counter+1};
        switch prop
            case 'zeroField'
                zeroField = val;
            case 'trialTypes'
                trialTypes = val;
            case 'outcomes'
                outcomes = val;                
            otherwise
        end
        counter=counter+2;
    end
    
%     trialTypes = unique(SessionData.TrialTypes);
%     outcomes = unique(SessionData.TrialOutcome);
    
    % initialize data fields
    data = struct(...
        'events', [],...
        'trials', []...
        );
    analysis.(event).data = repmat(data, max(trialTypes), max(outcomes));    
    
    
    for i = 1:length(trialTypes)
        type = trialTypes(i);
        for j = 1:length(outcomes)
            outcome = outcomes(j);
            trials = bpFilterTrials(SessionData, type, outcome);
            analysis.(event).data(type, outcome).trials = trials';
           % inialize events cell array
            events = cell(size(trials));
            analysis.(event).data(type, outcome).events=events;         
            for counter = 1:length(trials)
                trial = trials(counter);
                % always use first instance and start of zeroField for zeroing (each state has a start and end timestamp)                    
                if ~isempty(zeroField)
                    zeroTime = SessionData.RawEvents.Trial{trial}.States.(zeroField)(1,1);
                else
                    zeroTime = 0;
                end
                if isfield(SessionData.RawEvents.Trial{trial}.Events, event)
                    analysis.(event).data(type, outcome).events{counter} =...
                        SessionData.RawEvents.Trial{trial}.Events.(event) - zeroTime;
                end
            end
        end
    end
    