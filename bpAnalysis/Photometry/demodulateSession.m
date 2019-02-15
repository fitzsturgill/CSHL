function SessionData = demodulateSession(SessionData, varargin)
%% 1/17/2015 Fitz Sturgill
% demodulates the whole session and places it in nidaqData.demod field


%% default values
    defaults = {...
        'channels', [1 2];...
        'ACfilter', [0 0];...
        'refChannels', [1 2];...
        'lowpass', 15;...% corner frequency for lowpass filtering, default = 15Hz
        'forceAmp', 0;... % force demodulation even if the refChannel LED is off (i.e. it's amplitude = 0)
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

%% for very old sessions, sampleRate is hard coded (circa 2015)
try
    sampleRate = SessionData.Settings.sample_rate;
catch
    sampleRate = 6100;
end

% optionally emulate/apply AC Coupling to input data (in case photodiode was in DC coupling mode)
% first generate the transfer function coefficients
if isscalar(s.ACfilter)
    s.ACfilter = repmat(s.ACfilter, 1, length(s.channels));
end
if any(s.ACfilter)
    highCutoff = 15/(sampleRate / 2); % use a 30Hz filter cutoff, normalized relative to Nyquist frequency
    [z,p,k] = butter(5, highCutoff, 'high');
    [sos, g] = zp2sos(z,p,k);
end

% initialize demod cell array
SessionData.demod = cell(SessionData.nTrials, 2);



for counter = 1:length(s.channels)
    for trial = 1:SessionData.nTrials
        fCh = s.channels(counter);
%% Determine demodulation mode:    
%% case 1: reference data saved (prior to 6/2017)
        if trial == 1
            chIx = []; % figure out channel index below
            if ~isstruct(SessionData.NidaqData{1,2}) % figure this out from the first channel
                chIx = fCh; % for old style, you couldn't just acquire channel 2 (channel 1 was acquired too I believe, just didn't have useful data)
                if isfield(SessionData.TrialSettings(1,1), 'nidaq') % new style data acq parameters in nidaq field 
                    modF = SessionData.TrialSettings(1,1).nidaq.(['LED' num2str(fCh) '_f']); % this could change (I should store modF in settings not trialsettings)        
                else % old style
                    modF = SessionData.TrialSettings(1,1).(['LED' num2str(fCh) '_f']); % this could change (I should store modF in settings not trialsettings)        
                end
                demodMode = 1;
            else
%% case 2: reference data parameters saved as structure for each trial (introduced ~6/2017)
                refSettings = SessionData.NidaqData{trial, 2};
                % figure out channel index
                chIx = find(refSettings.channelsOn == fCh);
                modF = [];
                demodMode = 2;
                % find maximum LED amplitude on a given reference channel for force mode (where even if
                % LED is off you pretend it is on in order to determine the
                % level of bleedthrough from adjacent LEDs
                if s.forceAmp
                    amps = cellfun(@(x) x.amp(s.refChannels(counter)), SessionData.NidaqData(:,2));
                    forceAmp = max(amps);
                else
                    forceAmp = 0;
                end
            end            
        end
%% if data acq hiccupped and somehow didn't acquire during trial, replace with NaNs
        try
            assert(length(chIx) == 1); % channel should be present, otherwise throw error to fill with NaNs
            rawData = SessionData.NidaqData{trial,1}(:,chIx);
        catch 
            SessionData.demod{trial,fCh} = NaN(size(rawData)); % rawData from previous loop iteration
            SessionData.NidaqData{trial, 1}(:,fCh) = NaN(size(rawData));
            disp(['*** demodulateSession: Trial lacks NidaqData: # ' num2str(trial) ' ***']);
            continue
        end
        nSamples = size(rawData, 1);
        if s.ACfilter(counter) % whether or not to apply AC filtering for this channel
            rawData = filtfilt(sos, g, rawData);        
        end
        switch demodMode
            case 1
                refData = SessionData.NidaqData{trial,2}(1:nSamples,fCh);
                finalData = phDemod(rawData, refData, sampleRate, modF, s.lowpass); % lowpass corner freq                
            case 2
                refData = SessionData.NidaqData{trial,2};
                finalData = phDemod_v2(rawData, refData, s.refChannels(counter), sampleRate, 'forceAmp', forceAmp); 
        end
        SessionData.demod{trial,fCh} = finalData;
    end
end        


        