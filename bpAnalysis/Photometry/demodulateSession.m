function SessionData = demodulateSession(SessionData, lowpass, varargin)
%% 1/17/2015 Fitz Sturgill
% demodulates the whole session and places it in nidaqData.demod field

%%
if nargin < 2 || isempty(lowpass)
    lowpass = 15;  % corner frequency for lowpass filtering, default = 15Hz
end
%% default values
    defaults = {...
        'channels', [1 2];...
        'ACfilter', [0 0];...
        'refChannels', [1 2];...
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
if s.ACfilter
    highCutoff = 30/sampleRate * 2; % use a 30Hz filter cutoff
    [b, a] = butter(5, highCutoff, 'high');
end

% initialize demod cell array
SessionData.demod = cell(SessionData.nTrials, 2);


for trial = 1:SessionData.nTrials
    for counter = 1:length(s.channels)
        fCh = s.channels(counter);
%% Determine demodulation mode:    
%% case 1: reference data saved (prior to 6/2017)
        if trial == 1
            if ~isstruct(SessionData.NidaqData{1,2}) % figure this out from the first channel
                if isfield(SessionData.TrialSettings(1,1), 'nidaq') % new style data acq parameters in nidaq field 
                    modF = SessionData.TrialSettings(1,1).nidaq.(['LED' num2str(fCh) '_f']); % this could change (I should store modF in settings not trialsettings)        
                else % old style
                    modF = SessionData.TrialSettings(1,1).(['LED' num2str(fCh) '_f']); % this could change (I should store modF in settings not trialsettings)        
                end
                demodMode = 1;
            else
%% case 2: reference data parameters saved as structure for each trial (introduced ~6/2017)
                modF = [];
                demodMode = 2;
            end            
        end
%% if data acq hiccupped and somehow didn't acquire during trial, replace with NaNs
        try
            rawData = SessionData.NidaqData{trial,1}(:,fCh);
        catch 
            SessionData.demod{trial,fCh} = NaN(size(rawData)); % rawData from previous loop iteration
            SessionData.NidaqData{trial, 1}(:,fCh) = NaN(size(rawData));
            disp(['*** demodulateSession: Trial lacks NidaqData: # ' num2str(trial) ' ***']);
            continue
        end
        nSamples = size(rawData, 1);
        if s.ACfilter(counter) % whether or not to apply AC filtering for this channel
            rawData = filtfilt(b, a, rawData);        
        end
        switch demodMode
            case 1
                refData = SessionData.NidaqData{trial,2}(1:nSamples,fCh);
                finalData = phDemod(rawData, refData, sampleRate, modF, lowpass); % lowpass corner freq                
            case 2
                refData = SessionData.NidaqData{trial,2};
                finalData = phDemod_v2(rawData, refData, s.refChannels(counter), sampleRate); 
        end
        SessionData.demod{trial,fCh} = finalData;
    end
end        


        