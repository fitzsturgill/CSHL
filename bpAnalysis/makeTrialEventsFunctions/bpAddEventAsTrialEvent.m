function Event = bpAddEventAsTrialEvent(sessions, state, Event)

    if nargin < 3
        Event = {}; % going to build this across sessions by concatenation/ accretion
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
            trialStates = rawEvents{trial}.Events;
            if ~isfield(trialStates, state)
%                 break % break out of loop
%                 disp('**** WARNING in bpAddEventAsTrialEvent: are "missed" events (e.g. mouse did not lick ever during a trial) filled as NaNs or ommited for a trial s a field');
            else        
                theseStates{trial} = trialStates.(state);
            end
        end
        Event = [Event; theseStates]; % concatenate vertically
    end