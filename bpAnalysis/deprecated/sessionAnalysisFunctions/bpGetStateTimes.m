function ts = bpGetStateTimes(session, state, UniformOutput, start_vs_end)
    % state- name of state
%     UniformOutput- if 0, outputs cell array of size = nTrials, 1, default
%     = 1
%   start_vs_end-   valid for uniformoutupt = 1, first timestamp or last
    
    if nargin < 3
        UniformOutput = 1;
    end
    
    if nargin < 4
        start_vs_end = 'start';
    end
    
    if nargin == 4 && UniformOutput == 0
        error('!!! Error in bpGetStateTimes: first_vs_last is invalid when UniformOutput = 0 !!!');
    end
    
    if UniformOutput
        switch start_vs_end
            case 'start'
                ts = cellfun(@(x) x.States.(state)(1), session.SessionData.RawEvents.Trial, 'ErrorHandler', @trialLacksStateFcn)';
            case 'end'
                ts = cellfun(@(x) x.States.(state)(end), session.SessionData.RawEvents.Trial, 'ErrorHandler', @trialLacksStateFcn)';
        end
    else
        ts = cellfun(@(x) x.States.(state), session.SessionData.RawEvents.Trial, 'UniformOutput', 0, 'ErrorHandler', @trialLacksStateFcn)';
    end
    
    function sub = trialLacksStateFcn % this doesn't look to be necessary, I think that in RawEvents.Trial, non-visited events are filled with NaNs
        if UniformOutput
            sub = NaN; % substitute NaNs
        else
            sub = {}; % substitute {};
        end
    end
        
end

        