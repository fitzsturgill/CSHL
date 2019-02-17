function sessions = demodulateSessions(sessions, lowpass, varargin)
    
    if nargin < 2
        lowpass = []; % use default from demodulateSession
    end
    
    if nargin < 3
        varargin = {}; % again use defaults
    end

    for counter = 1:length(sessions)
        disp(['*** Demodulating: ' sessions(counter).filename '***']);
        sessions(counter).SessionData = demodulateSession(sessions(counter).SessionData, lowpass, varargin);
    end
    