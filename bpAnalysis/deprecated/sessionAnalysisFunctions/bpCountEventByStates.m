function eventCount = bpCountEventByStates(session, event, referenceState, varargin)
    % if end
    defaults = {...
        'event', event;...        
        'referenceState', referenceState;... % event count is started at beginning of referenceState by default
        'referenceFromEnd', 0;...
        'window', [0 0];... % with respect to referenceState, limit or extend the window defined by state(s) boundaries       
%        window(1) offsets start of count from         
        'endState', [];... % if endState is provided, stop count at window(2) or end of endState, whichever comes first.
        };
    
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    if s.referenceFromEnd && isempty(s.endState) && isempty(find(s.window))
        error('Error in bpCountEventByStates, provide endState or ending window');
    end
    nTrials = session.SessionData.nTrials;
    RawEvents = session.SessionData.RawEvents;
    eventCount = struct(...
        'count', zeros(1, nTrials),...
        'duration', zeros(1, nTrials),...
        'rate', zeros(1, nTrials),...
        'settings', repmat(s, 1, nTrials)...
        );
    
    
    for trial = 1:nTrials
%         Get start and stop time
        if isfield(RawEvents.Trial{trial}.States, s.referenceState)
            % determine t1 (start time)
            if ~s.referenceFromEnd
                t1 = RawEvents.Trial{trial}.States.(s.referenceState)(1) + s.window(1); % remember window(1) can be negative!, this is useful
            else
                t1 = RawEvents.Trial{trial}.States.(s.referenceState)(end) + s.window(1); % remember window(1) can be negative!, this is useful                
            end
            % determine t2 (end time)
            if ~isempty(s.endState)  && isfield(RawEvents.Trial{trial}.States, s.endState)
                if s.window(2) == 0 % if terminating window not provided default to end of endState
                    t2 = RawEvents.Trial{trial}.States.(s.endState)(end);
                else % if end window provided, use whichever comes first (window end or endState end)
                    t2 = min(t1 + s.window(2), RawEvents.Trial{trial}.States.(s.endState)(end));
                end
            elseif  s.window(2) == 0 % endState not provided, end window not provided
                t2 = RawEvents.Trial{trial}.States.(s.referenceState)(end); % use end of referenceState
            else % endState not provided, end window provided
                t2 = t1 + s.window(2); % window(2) offset from start of referenceState
            end
            if t2 < t1
                warning('Invalid window in bpCountEventByStates'); % count should still be 0
                duration = NaN;
            else
                duration = t2 - t1;
            end
            
            if isfield(RawEvents.Trial{trial}.Events, s.event)
                te = RawEvents.Trial{trial}.Events.(s.event);
                count = sum(te >= t1 & te <= t2);
            else
                count = 0;
            end
            rate = count/duration;
        else
            count = NaN;
            duration = NaN;
            rate = NaN;
        end
        eventCount.count(trial) = count;
        eventCount.rate(trial) = rate;
        eventCount.duration(trial) = duration;
    end
%% OLD
% function analysis = bpCountEventByStates(analysis, event, label, states, varargin)
%     if ~isfield(analysis, event)
%         disp(['*** Error in countEventByStates: Analysis for event: ' event ' does not yet exist ***']);
%         return
%     end
%     
%     if ischar(states)
%         states = {states};
%     end
%     
%     %% default values
%     trialTypes = unique(SessionData.TrialTypes);
%     outcomes = unique(SessionData.TrialOutcome);    
%     
%     % as a default, try and supply zeroField from settings
% %     e.g. for session.analysis.Port1In.settings, ans = {'zeroField', 'DeliverStimulus'}
%     zfi = find(cellfun(@(x) strcmp(x, 'zeroField'), analysis.(event).settings)); % zerofield index
%     if zfi
%         zeroField = analysis.(event).settings{zfi + 1};
%     else
%         zeroField = '';
%     end
% 
%     %parse input parameter pairs
%     counter = 1;
%     while counter+1 <= length(varargin) 
%         prop = varargin{counter};
%         val = varargin{counter+1};
%         switch prop
%             case 'zeroField'
%                 zeroField = val;
%             case 'trialTypes'
%                 trialTypes = val;
%             case 'outcomes'
%                 outcomes = val;                
%             otherwise
%         end
%         counter=counter+2;
%     end
% %% if ~isfield(analysis, event)
%         disp(['*** Error in countEventByStates: Analysis for event: ' event ' does not yet exist ***']);
%         return
%     end
%     
%     if ischar(states)
%         states = {states};
%     end
%     
%     %% default values
%     trialTypes = unique(SessionData.TrialTypes);
%     outcomes = unique(SessionData.TrialOutcome);    
%     
%     % as a default, try and supply zeroField from settings
% %     e.g. for session.analysis.Port1In.settings, ans = {'zeroField', 'DeliverStimulus'}
%     zfi = find(cellfun(@(x) strcmp(x, 'zeroField'), analysis.(event).settings)); % zerofield index
%     if zfi
%         zeroField = analysis.(event).settings{zfi + 1};
%     else
%         zeroField = '';
%     end
% 
%     %parse input parameter pairs
%     counter = 1;
%     while counter+1 <= length(varargin) 
%         prop = varargin{counter};
%         val = varargin{counter+1};
%         switch prop
%             case 'zeroField'
%                 zeroField = val;
%             case 'trialTypes'
%                 trialTypes = val;
%             case 'outcomes'
%                 outcomes = val;                
%             otherwise
%         end
%         counter=counter+2;
%     end
% %%

