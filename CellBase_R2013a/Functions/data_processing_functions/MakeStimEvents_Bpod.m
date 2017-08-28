function MakeStimEvents_Bpod(sessionpath,varargin)

%% FS MOD of Balazs' makeStimEvents function, relies solely on Neuralynx TTLs to capture stim events
%MAKESTIMEVENTS2   Create stimulus events structure.
%   MAKESTIMEVENTS2(SESSIONPATH) constructs and saves stimulus events
%   structure ('StimEvents') for the stimulation session in the folder
%   SESSIONPATH. Stimulation pulses are detected as TTLs recorded by the
%   data acquisition system (e.g. Neuralynx) and their number is verified
%   by the log file of the stimulation control software. Note however, that
%   the recorded TTLs are exclusive source of stimulus information. Note
%   also, that the Neuralynx Event file should first be converted to Matlab
%   .mat file before MAKESTIMEVENTS2 can run.
%
%   Stimulus time stamps (e.g. pulse onset and offset, burst onset and
%   offset, protocol start and end), time intervals (e.g. inter-burst
%   interval, burst duration, pulse duration), stimulation-related
%   statistics (e.g. frequency) are extracted/calculated from the recorded
%   TTL time stamps and event strings. The data are loaded in a structure
%   and saved under the name 'StimEvents.mat'.
%
%   By using the MAKESTIMEVENTS2(SESSIONPATH,PARAMETER1,VALUE1,...) syntax, 
%   the TTLs actually used during the recording as well as the name of
%   stimulus protocol file can be specified (default values provided):

%       'PulseNttl', 16384 - TTL signalling pulse onset and offset time
%           (TTL split between stimulation device and data acquisition
%           system)
%
%   See also PARSETTLS.

%   Edit log: BH 4/2/13

% Parse input arguments
% default_args={...
%     'PulseNttl',           16384;... % TTL signalling pulse onset and offset time
%     'PulsePort', 
%     };
% [g, error] = parse_args(default_args,varargin{:});

prs = inputParser;
addRequired(prs,'sessionpath',@ischar)  % pathname for session
addParamValue(prs,'PulseNttl',128,@isnumeric)   % TTL signalling pulse onset and offset time
addParamValue(prs, 'PulsePort', 0, @isnumeric) % Port number of PulseNttl
% addParamValue(prs, 'BpodPort', 1, @isnumeric)
% addParamValue(prs,'BurstStartNttl',64,@isnumeric)   % TTL for start of burst stimulation
% addParamValue(prs,'ProtocolStartNttl',128,@isnumeric)   % TTL sent at the beginning of each stimulation protocol

%           of simulation and behavior protocols (cell array); no action is
%           implemented in the current version for behavior protocols 
parse(prs,sessionpath,varargin{:})
g = prs.Results;

% Check session path
if ~isdir(sessionpath)
    error('Session path is not valid.');
end
cd(sessionpath)


% Load Neuralynx events
try
    load('EVENTS.mat')
catch
    error('EVENTS file not found; make sure Nlx files have been converted.');
end





% % Pulse TTL's
% [parsed_ttls onttl offttl] = parseTTLs(Events_Nttls);   %#ok<ASGLU> % parse TTL's
% inx = size(onttl,2) - log2(g.PulseNttl);
% pepon = find(onttl(:,inx));
% pepoff = find(offttl(:,inx));
% Epon = intersect(find(Events_EventIDs==g.PulseEventID),pepon);   % pulse on
% Epoff = intersect(find(Events_EventIDs==g.PulseEventID),pepoff);   % pulse off
%% FS MOD
% if ~isempty(g.PulsePort)
    PortID = eventPortFromEventStrings(Events_EventStrings);
    Pulses = PortID == g.PulsePort & Events_Nttls == g.PulseNttl;
    TrialStart = PortID == g.PulsePort & Events_Nttls == g.PulseNttl;
%     Epon = intersect(find(PortID == g.PulsePort), Epon);
%     Epoff = intersect(find(PortID == g.PulsePort), Epoff);   
% end
% BSInfo = find(Events_EventIDs==3&Events_Nttls==g.BurstStartNttl);  % burst start info indices

% %% Count bursts and pulses for each protocol
% [valid_bs valid_pulses] = deal(cell(NumProt,1));
% [num_bs num_pulses] = deal(zeros(NumProt,1));
% for iProt = 1:NumProt   % loop through the protocols
%     pstart = Events_TimeStamps(ProtEvents(1,iProt));   % start time
%     pend = Events_TimeStamps(ProtEvents(2,iProt));   % end time
%     valid_bs{iProt} = BSInfo(Events_TimeStamps(BSInfo)>pstart&Events_TimeStamps(BSInfo)<pend);   % burst onset indices for the protocol
%     valid_pulses{iProt} = Epon(Events_TimeStamps(Epon)>pstart&Events_TimeStamps(Epon)<pend);    % pulse onset indices for the protocol
%     num_bs(iProt) = length(valid_bs{iProt});   % number of bursts
%     num_pulses(iProt) = length(valid_pulses{iProt});   % number of pulses
%     try   % cross-check number of bursts in the saved protocol file
%         SP = load([sessionpath filesep PName{iProt}]);   % load protocol file
%         if ~length(find(SP.TrialsCompleted))==length(valid_bs)
%             warning('MakeStimEvents:burstNoMismatch','Number of recorded burst events does not match the protocol file.');
%         end
%     catch
%         warning('MakeStimEvents:loadFailure','Unable to load protocol file ''%s''.',PName{iProt});
%     end
% end
% nvalid_bs = sum(num_bs(logical(val_stim_protocols)));  % total number of bursts for all stim. protocols
% nvalid_pulses = sum(num_pulses(logical(val_stim_protocols)));   % total number of pulses for all stim. protocols

% Preallocate stimulus events
nvalid_pulses = sum(Pulses);
SE = struct(...
    'TrialStart', zeros(1, nvalid_pulses),...
    'Pulse', nan(1,nvalid_pulses)...
    );

SE.Pulse = Events_TimeStamps(Pulses);

% SE = struct;
% SE.ProtocolName = cell(1,nvalid_pulses); % e.g. LaserStimProtocol2
% SE.ProtocolID = cell(1,nvalid_pulses);   % 'S' for stim. and 'B' for behav. protocols
% SE.ProtocolStart = nan(1,nvalid_pulses); % time stamp for protocol begin
% SE.ProtocolEnd = nan(1,nvalid_pulses);   % time stamp for protocol end
% SE.StimType = nan(1,nvalid_pulses);      % 1 for single pulse, 2 for burst stimualtion
% SE.TrialsCompleted = nan(1,nvalid_pulses); % not implemented
% SE.TrialStart = zeros(1,nvalid_pulses);  % not implemented
% 
% SE.BurstOn = nan(1,nvalid_pulses);       % onset of first pulse in the burst
% SE.BurstOff = nan(1,nvalid_pulses);      % offset of last pulse in the burst
% SE.BurstDur = nan(1,nvalid_pulses);      % duration of burst (between last pulse offset and first pulse onset
% SE.BurstIBI = nan(1,nvalid_pulses);      % inter-burst interval
% SE.BurstNPulse = nan(1,nvalid_pulses);   % number of pulses in burst (Events_EventStrings)
% SE.PrevBurstOff = nan(1,nvalid_pulses);  % end of previous burst
% SE.NextBurstOn = nan(1,nvalid_pulses);   % start of next burst
% SE.PreBurstIBI = nan(1,nvalid_pulses);   % time since PrevBurstOff
% SE.PostBurstIBI = nan(1,nvalid_pulses);  % time to NextBurstOn
% SE.BurstID = nan(1,nvalid_pulses);       % burst rank in the protocol (NTrial)
% 
% SE.PulseOn = nan(1,nvalid_pulses);       % onset of stim. pulse
% SE.PulseOff = nan(1,nvalid_pulses);      % offset of stim. pulse
% SE.PulseDur = nan(1,nvalid_pulses);      % duration of light pulse (Events_EventStrings)
% SE.PulseIPI = nan(1,nvalid_pulses);      % inter-pulse interval (Events_EventStrings)
% SE.PulseFreq = nan(1,nvalid_pulses);     % pulse frequency (Events_EventStrings)
% SE.PulsePower = nan(1,nvalid_pulses);    % stimulus intensity (Events_EventStrings)
% SE.PrevPulseOff = nan(1,nvalid_pulses);  % offset of previous pulse
% SE.NextPulseOn = nan(1,nvalid_pulses);   % onset of next pulse
% SE.PrePulseIPI = nan(1,nvalid_pulses);   % time since PrevPulseOff
% SE.PostPulseIPI = nan(1,nvalid_pulses);  % time to NextPulseOn
% 
% SE.FirstPulse = nan(1,nvalid_pulses);    % is this the first pulse in a burst?
% SE.ZeroPulse = nan(1,nvalid_pulses);     % extrapolate one pulse time back
% SE.LastPulse = nan(1,nvalid_pulses);     % is this the last pulse in a burst?
% SE.PulseNum = nan(1,nvalid_pulses);      % pulse rank in the burst
% 
% % Fill stimulus events structure
% pco = 1;  % pulse counter
% for iProt = 1:NumProt   % loop through the protocols
%     if val_stim_protocols(iProt) == 1   % if it is a stimulation protocol
%         bursts = valid_bs{iProt};   % burst indices
%         pulses = valid_pulses{iProt};    % pulse indices
%         for iB = 1:length(bursts)   % loop through bursts
%             if iB == length(bursts)   % last burst - use protocol end time stamp
%                 pulses_ind = find(Events_TimeStamps(pulses)>Events_TimeStamps(bursts(iB))...
%                     &Events_TimeStamps(pulses)<Events_TimeStamps(ProtEvents(2,iProt)));  % find pulse indices that are within the burst
%             else
%                 pulses_ind = find(Events_TimeStamps(pulses)>Events_TimeStamps(bursts(iB))...
%                     &Events_TimeStamps(pulses)<Events_TimeStamps(bursts(iB+1))); % find pulse indices that are within the burst
%             end
%             eval(Events_EventStrings{bursts(iB)});
%             if ~exist('BurstPulseIPI','var')
%                 BurstDur = 2;   % burst duration not saved in early protocol files; it was kept 2s
%                 disp('2 sec burst duration assumed.')
%                 BurstPulseIPI = BurstDur / BurstNPulse - BurstPulseDur;  % burst pulse IPI not saved in early protocols (2s burst duration assumed) - will become obsolete
%             end
%             if BurstNPulse ~= length(pulses_ind),
%                 warning('MakeStimEvents:pulseNoMismatch','No. of pulses does not match between recording and behavior systems for burst no.%d and ProtocolID %s',iB,PName{iProt});
%             end
%             for iP = 1:length(pulses_ind)   % loop through pulses
%                 if length(pulses_ind) == 1   % stimulation type (single pulse/burst)
%                     SE.StimType(pco) = 1;
%                 else
%                     SE.StimType(pco) = 2;
%                 end
%                 SE.ProtocolName{pco} = PName{iProt};   % protocol name
%                 SE.ProtocolID{pco} = PID{iProt};   % protocol ID
%                 if exist('NTrial','var')  % burst rank (serial no.)
%                     SE.BurstID(pco) = NTrial;
%                 else
%                     SE.BurstID(pco) = iB;   % use cycle variable if the info was not saved in the recording system
%                 end
%                 SE.PulseOn(pco) = Events_TimeStamps(pulses(pulses_ind(iP)));   % pulse onset: from recorded pulse TTL
%                 SE.PulseOff(pco) = Events_TimeStamps(Epoff(find(Epoff>pulses(pulses_ind(iP)),1,'first')));   % pulse off: first pulse offset after pulse onset
%                 SE.PulseDur(pco) = BurstPulseDur; % duration of stim. pulse (Events_EventStrings)
%                 SE.PulseIPI(pco) = BurstPulseIPI; % inter-pulse interval (Events_EventStrings)
%                 SE.PulseFreq(pco) = 1 / (BurstPulseDur + BurstPulseIPI); % pulse frequency
%                 SE.PulsePower(pco) = BurstPulsePower; % stim. intensity (Events_EventStrings)
%                 SE.BurstNPulse(pco) = BurstNPulse; % number of pulses in burst (Events_EventStrings)
%                 SE.PulseNum(pco) = iP;  % pulse rank in the burst
%                 if iP == 1   % first pulse
%                     SE.FirstPulse(pco) = 1;
%                     SE.ZeroPulse(pco) = SE.PulseOn(pco) - (BurstPulseDur + BurstPulseIPI);  % extrapolate one pulse time back
%                     SE.BurstOn(pco) = SE.PulseOn(pco);   % burst onset
%                     if iB==1
%                         SE.ProtocolStart(pco) = Events_TimeStamps(PS(iProt));  % protocol start
%                     end 
%                 end
%                 if iP == length(pulses_ind)   % last pulse
%                     SE.LastPulse(pco) = 1;
%                     SE.BurstOff(pco) = SE.PulseOff(pco);  % burst end
%                     if iB == length(bursts)   % last burst
%                         SE.ProtocolEnd(pco) = SE.PulseOff(pco);  % protocol end: last pulse offset
%                     end
%                 end
%                 
%                 pco = pco+1;   % update pulse counter
%             end
%         end
%     else
%         % no action implemented for behavior protocols
%     end
% end
% 
% SE.PrevPulseOff = [SE.ProtocolStart(1) SE.PulseOff(1:end-1)];   % previous pulse offset (first protocol start for first pulse)
% SE.NextPulseOn = [SE.PulseOn(2:end) SE.ProtocolEnd(end)];   % next pulse onset (last protocol end for last pulse)
% SE.PrevBurstOff = [SE.ProtocolStart(1) SE.BurstOff(1:end-1)];   % previous burst offset (first protocol start for first pulse)
% SE.NextBurstOn = [SE.BurstOn(2:end) SE.ProtocolEnd(end)];   % next burst onset (last protocol end for last pulse)
% SE.PrePulseIPI = SE.PulseOn - SE.PrevPulseOff;   % inter-pulse interval before pulse
% SE.PostPulseIPI = SE.PulseOff - SE.NextPulseOn;   % inter-pulse interval after pulse
% SE.PreBurstIBI = SE.BurstOn - SE.PrevBurstOff;   % previous inter-burst interval
% SE.PostBurstIBI = SE.BurstOff - SE.NextBurstOn;   % next inter-burst interval
% SE.BurstDur = SE.BurstOff - SE.BurstOn;   % burst duration (real, not from EventStrings)
% 
% % Addition for 'omit protocol'
% freqs = unique(SE.PulseFreq);  % different stim. frequencies
% omitpulse = nan(size(SE.PulseOn));
% for k = 1:length(freqs)   % loop through frequencies
%     vo = 2 / freqs(k);     % length of ISI if omission occurs
%     dfp = diff(SE.PulseOn);
%     prec = 0.003;  % precision
%     omitinx = find(dfp>vo-prec&dfp<vo+prec&SE.PulseFreq(1:end-1)==freqs(k));
%     omitpulse(omitinx) = SE.PulseOn(omitinx) + vo / 2;   % time point of the omitted pulse
% end
% SE.OmitPulse = omitpulse;  % time stamps for omitted pulses

% Save 'StimEvents'
save([sessionpath filesep 'StimEvents.mat'],'-struct','SE')