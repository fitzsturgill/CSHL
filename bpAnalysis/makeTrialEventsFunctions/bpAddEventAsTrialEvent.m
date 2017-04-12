function TrialEvent = bpAddEventAsTrialEvent(sessions, eventName, TrialEvent)

    if nargin < 3
        TrialEvent = {}; % going to build this across sessions by concatenation/ accretion
    end
    % eventName- string corresponding to name of a given Bpod event
    % TrialEvent- optional, if you want to add onto previouly generated
    % stateEvent cell array
    for si = 1:length(sessions)
        session = sessions(si);
        nTrials = session.SessionData.nTrials;
        rawEvents = session.SessionData.RawEvents.Trial; % cell array
        theseEvents = repmat({NaN}, nTrials, 1); % initialize for each session, NaNs by default
        for trial = 1:nTrials
            trialStates = rawEvents{trial}.Events;
            if ~isfield(trialStates, eventName)
            else        
                theseEvents{trial} = trialStates.(eventName);
            end
        end
        TrialEvent = [TrialEvent; theseEvents]; % concatenate vertically
    end