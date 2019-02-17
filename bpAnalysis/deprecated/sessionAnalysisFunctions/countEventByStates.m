function analysis = countEventByStates(analysis, event, epochLabel, states, varargin)
    
    if ~isfield(analysis, event)
        disp(['*** Error in countEventByStates: Analysis for event: ' event ' does not yet exist ***']);
        return
    end
    
    if ischar(states)
        states = {states};
    end
    
    %% default values
    trialTypes = unique(SessionData.TrialTypes);
    outcomes = unique(SessionData.TrialOutcome);    
    
    % as a default, try and supply zeroField from settings
%     e.g. for session.analysis.Port1In.settings, ans = {'zeroField', 'DeliverStimulus'}
    zfi = find(cellfun(@(x) strcmp(x, 'zeroField'), analysis.(event).settings)); % zerofield index
    if zfi
        zeroField = analysis.(event).settings{zfi + 1};
    else
        zeroField = '';
    end

    %parse input parameter pairs
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
%%

