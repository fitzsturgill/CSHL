function eventCount = bpCountEventByStates2(sessions, event, referenceState, varargin)
    % if end
    defaults = {...
        'event', event;...    
        'eventCount', [];... % create new eventCount by default (don't append)
        'referenceState', referenceState;... % event count is started at beginning of referenceState by default
        'referenceFromEnd', 0;...
        'window', [0 0];... % with respect to referenceState, limit or extend the window defined by state(s) boundaries       
%        window(1) offsets start of count from         
        'endState', [];... % if endState is provided, stop count at window(2) or end of endState, whichever comes first.
        };
    
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    if s.referenceFromEnd && isempty(s.endState) && isempty(find(s.window))
        error('Error in bpCountEventByStates, provide endState or ending window when starting at end of reference state');
    end
    
    % if eventCount provided, append otherwise, create new eventCount
    % structure
    if isempty(s.eventCount)        
        eventCount = struct(...
            'count', [],...  % easier to just concatenate eventCount across sessions (build by accretion) rather than initializing in beginning
            'duration', [],...
            'rate', [],...
            'settings', []...
            );
    else
        eventCount = s.eventCount;
    end
    
    for si = 1:length(sessions) % si = session index
        session = sessions(si);
        nTrials = session.SessionData.nTrials;
        RawEvents = session.SessionData.RawEvents;
        count = zeros(nTrials, 1);
        duration = zeros(nTrials, 1);
        rate = zeros(nTrials, 1);
        settings = repmat(s, nTrials, 1);
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
                    tduration = NaN;
                else
                    tduration = t2 - t1;
                end

                if isfield(RawEvents.Trial{trial}.Events, s.event)
                    te = RawEvents.Trial{trial}.Events.(s.event);
                    tcount = sum(te >= t1 & te <= t2);
                else
                    tcount = 0;
                end
                trate = tcount/tduration;
            else
                tcount = NaN;
                tduration = NaN;
                trate = NaN;
            end
            count(trial) = tcount;
            rate(trial) = trate;
            duration(trial) = tduration;
        end
        eventCount.count = [eventCount.count; count];
        eventCount.rate = [eventCount.rate; rate];
        eventCount.duration = [eventCount.duration; duration];
        eventCount.settings = [eventCount.settings; settings];
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

