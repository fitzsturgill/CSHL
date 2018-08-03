function trialStartTimes = getBehaviorStartTimes(Nttls, EventStrings, TimeStamps, varargin)
    % select trial start times that don't include a laser pulse    
    % 8/2018, modified to include ALL Bpod trial starts such that correct
    % trial starts are identified by time shift method
    
        %% optional parameters, first set defaults
    defaults = {...
        'laserPort', 0;... % which neuralynx port laser is plugged into
        'laserNttl', 128;... % laser Nttl
        'behaviorPort', 1;... % bpod sync port connected
        'behaviorNttl', 128;... % bpod trial start Nttl
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    s.laserPort = 0;
    s.laserNttl = 128;
    s.behaviorPort = 1;
    s.behaviorNttl = 128;
    PortID = eventPortFromEventStrings(EventStrings);
    allStartsIx = PortID == s.behaviorPort & Nttls == s.laserNttl;
    allStarts = TimeStamps(allStartsIx);
%     laserIx = PortID == s.laserPort & Nttls == s.laserNttl;
%     pulses = TimeStamps(laserIx);
    trialStartTimes = allStarts;
    
    
%     behaviorStartsIx = false(size(allStarts));

%     for counter = 1:length(allStarts)
%         thisStart = allStarts(counter);
%         
% 
%         if counter == length(allStarts)
%             nextStart = Inf;
%         else
%             nextStart = allStarts(counter + 1);
%         end
%         if any((pulses > thisStart) & (pulses < nextStart)) % if there are pulses in between a trial and the next, it's a laser tagging trial so skip
%             continue
%         end
%         behaviorStartsIx(counter) = true;        
%     end
%     
%     trialStartTimes = allStarts(behaviorStartsIx);
    