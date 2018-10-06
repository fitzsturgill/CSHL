function sessions = bpLoadSessionsFromTE(TE)
    % load session data from which TE is derived
        if isfield(TE, 'sessions')
            sessions = bpLoadSessions([], {TE.sessions.filename}, {TE.sessions.filepath});
        else
            error('*** TE doesn''t contain sessions field ***');
        end
    
    