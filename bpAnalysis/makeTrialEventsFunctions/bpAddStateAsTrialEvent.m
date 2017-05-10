function stateEvent = bpAddStateAsTrialEvent(sessions, state, stateEvent, matchMode)
    % state- string corresponding to name of a given Bpod state, trailing *
    % sets matchMode == 1
    % stateEvent- optional, if you want to add onto previuosly generated
    % stateEvent cell array    
    % matchMode, if set to 1 it groups states matching state string (in
    % case you don't want to add a trailing * or want to group states
    % containing rather than starting with 'state'

    if nargin < 3 || isempty(stateEvent)
        stateEvent = {}; % going to build this across sessions by concatenation/ accretion
    end
    
    if nargin < 4
        matchMode = 0;
    end

    
    if state(end) == '*'
        matchMode = true;
        state = state(1:end-1);
    end
    
    for si = 1:length(sessions)
        session = sessions(si);
        nTrials = session.SessionData.nTrials;
        rawEvents = session.SessionData.RawEvents.Trial; % cell array
        theseStates = repmat({NaN}, nTrials, 1); % initialize for each session, NaNs by default
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
            theseStates{trial} = timeStamps;
        end
        stateEvent = [stateEvent; theseStates]; % concatenate vertically
    end