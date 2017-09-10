function TE = makeTE_CuedOutcome_Odor_Complete_Nlx

    % add automatic loading later
    sessions = bpLoadSessions; % select Bpod session file from within CellBase file structure, e.g. within C:\FitzData\Cellbase_dev\CD4\170909a
    
    if length(sessions) > 1
        error('*** only handles single sessions for now ***');
    end
    TE = makeTE_CuedOutcome_Odor_Complete(sessions);
    
    load(fullfile(sessions.filepath, 'Events.mat'));   % load Neuralynx events
    % now Events_EventIDs, Events_EventStrings, Events_Extras,
    % Events_Nttls, and Events_TimeStamps are loaded into the local
    % scope/workspace
    
    