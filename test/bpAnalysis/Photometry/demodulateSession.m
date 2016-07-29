function SessionData = demodulateSession(SessionData, lowpass, varargin)
    %% 1/17/2015 Fitz Sturgill
    % demodulates the whole session and places it in nidaqData.demod field
    
    %%
    if nargin < 2 || isempty(lowpass)
        lowpass = 15;  % corner frequency for lowpass filtering, default = 15Hz
    end
    %% default values
    ACfilter = 0; 
    
    % parse input parameter pairs
    counter = 1;
    while counter+1 <= length(varargin) 
        prop = varargin{counter};
        val = varargin{counter+1};
        switch prop
            case 'ACfilter'
                ACfilter = val;
            otherwise
        end
        counter=counter+2;
    end    
    %%
    
    try
        sampleRate = SessionData.Settings.sample_rate;
    catch
        sampleRate = 6100;
    end
    nChannels = 2;
    
    % optionally emulate/apply AC Coupling to input data (in case photodiode was in DC coupling mode)
    % first generate the transfer function coefficients
    if ACfilter
        highCutoff = 30/sampleRate * 2; % use a 30Hz filter cutoff
        [b, a] = butter(5, highCutoff, 'high');
    end
    
    % initialize demod cell array
    SessionData.demod = cell(SessionData.nTrials, 2);
    
    for fCh = 1:nChannels
        if isfield(SessionData.TrialSettings(1,1), 'nidaq') % new style data acq parameters in nidaq field 
            modF = SessionData.TrialSettings(1,1).nidaq.(['LED' num2str(fCh) '_f']); % this could change (I should store modF in settings not trialsettings)        
        else % old style
            modF = SessionData.TrialSettings(1,1).(['LED' num2str(fCh) '_f']); % this could change (I should store modF in settings not trialsettings)        
        end
        for trial = 1:SessionData.nTrials
            try
                rawData = SessionData.NidaqData{trial,1}(:,fCh);
            catch % if data acq hiccupped and somehow didn't acquire during trial
%                 attempty to replace data with NaNs of size from previous
%                 trial
                SessionData.demod{trial,fCh} = NaN(size(rawData, 1), nChannels); % rawData from previous loop iteration
                SessionData.NidaqData{trial, 1}(:,fCh) = NaN(size(rawData));
                disp(['*** demodulateSession: Trial lacks NidaqData: # ' num2str(trial) ' ***']);
            end
            nSamples = size(rawData, 1);
            refData = SessionData.NidaqData{trial,2}(1:nSamples,fCh);
    %         [finalData, dmf, dmp] = decode_lockin_fn(rawData, refData, [], modF, sampleRate, 0);
            if ACfilter
                rawData = filtfilt(b, a, rawData);        
            end 
            finalData = phDemod(rawData, refData, sampleRate, modF, lowpass); % lowpass corner freq
            SessionData.demod{trial,fCh} = finalData;
        end
    end        
        
        
        