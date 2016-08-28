function stateEvent = bpAddStateAsTrialEvent(sessions, state, stateEvent)

    if nargin < 3
        stateEvent = {}; % going to build this across sessions by concatenation/ accretion
    end
    % stateEvent- string corresponding to name of a given Bpod state
    % states- optional, if you want to add onto previuosly generated
    % stateEvent cell array
    for si = 1:length(sessions)
        session = sessions(si);
        nTrials = session.SessionData.nTrials;
        rawEvents = session.SessionData.RawEvents.Trial; % cell array
        theseStates = repmat({NaN}, nTrials, 1); % initialize for each session, NaNs by default
        for trial = 1:nTrials
            trialStates = rawEvents{trial}.States;
            if ~isfield(trialStates, state)
                break % break out of loop
            else        
                theseStates{trial} = trialStates.(state);
            end
        end
        stateEvent = [stateEvent; theseStates]; % concatenate vertically
    end