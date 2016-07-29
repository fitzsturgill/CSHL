function SessionData = bpExcludeTrials(SessionData, exclude, bool, reset)

    if nargin < 3
        bool = 0; % default is to exclude trials, optionally set to 1 to include
    end
    
    if nargin < 4
        reset = 0; %  
    end

    if ~isfield(SessionData, 'IncludeTrials') || reset
        SessionData.IncludeTrials = ones(1, SessionData.nTrials);
    end
%     exclude = exclude(1:min(max(exclude), SessionData.nTrials)); % make sure you haven't accidently tried to exclude trials that don't exist
    SessionData.IncludeTrials(exclude) = bool;
    
    