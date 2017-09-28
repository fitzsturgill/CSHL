function stateEvent = bpAddStateAsTrialEvent(sessions, state, stateEvent, matchMode, whichTimeStamps)
    % state- string corresponding to name of a given Bpod state, trailing *
    % sets matchMode == 1
    % stateEvent- optional, if you want to add onto previuosly generated
    % stateEvent cell array    
    % matchMode, if set to 1 it groups states matching state string (in
    % case you don't want to add a trailing * or want to group states
    % containing rather than starting with 'state' 
    % whichTimeStamp (New Feature 9/2017): useful for cellbase because
    % event triggers per trial are normally scalar or multiple independent
    % events (e.g. lick times)
    % whichTimeStamp (New Feature 9/2017)- default 'all', other options: 'first' 'last'
    % if 'all', creates a nTrials x 1 cell array containing all the time stamps for
    % that state in a given trial. Contents of each element correspond to a
    % nTimesStateVisited x 2 matrix, with columns containing the start end
    % end values for each state visitation within a trial
    % 'first'- creates a nTrials x 1 double filled with the first time
    % stamp
    % 'last'-  creates a nTrials x 1 double filled with the last time stamp
    % for a given trial


    if nargin < 5
        whichTimeStamps = 'all'; 
    end
    
    if nargin < 3 || isempty(stateEvent)
        if strcmp(whichTimeStamps, 'all')
            stateEvent = {}; % going to build this across sessions by concatenation/ accretion
        else
            stateEvent = [];
        end
    end
    
    if nargin < 4
        matchMode = 0;
    end
    
    % cho
    if nargin < 5
        whichTimeStamps = 'all'; 
    end

    
    if state(end) == '*'
        matchMode = true;
        state = state(1:end-1);
    end
    
    for si = 1:length(sessions)
        session = sessions(si);
        nTrials = session.SessionData.nTrials;
        rawEvents = session.SessionData.RawEvents.Trial; % cell array
        if strcmp(whichTimeStamps, 'all')
            theseStates = repmat({NaN}, nTrials, 1); % initialize for each session, NaNs by default
        else
            theseStates = NaN(nTrials, 1);
        end
        for trial = 1:nTrials
            trialStates = rawEvents{trial}.States;
%             strncmp('Reward', fitz2, length('Reward'))
            if matchMode
                timeStamps = [];
                trialFields = fieldnames(trialStates);
                for counter = 1:length(trialFields)
                    field = trialFields{counter};
                    if strfind(field, state)
                        timeStamps = [timeStamps; trialStates.(field)];
                    end
                end
                if isempty(timeStamps)
                    continue % to retain NaN
                else
                    timeStamps = sortrows(timeStamps);
                end
            else
                if ~isfield(trialStates, state)
                    continue % to retain NaN
                else        
                    timeStamps = trialStates.(state);
                end
            end
            switch whichTimeStamps
                case 'all'
                    theseStates{trial} = timeStamps;
                case 'first'
                    theseStates(trial) = timeStamps(1);
                case 'last'
                    theseStates(trial) = timeStamps(end);
                otherwise
            end            

        end
        stateEvent = [stateEvent; theseStates]; % concatenate sessions vertically
    end