function exportPupil
    % export pupilometry to .mat file
    global state
    
    
    
    
    Photometry = struct(...
        'data', [],...          % length = number of Channels
        'settings', s,...
        'sampleRate', sampleRate,... % downsampled sample rate
        'startTime',startTimes...
    );
    
    data = struct(...
        'dFF', [],... %deltaF/F
        'raw', [],... %
        'blF', [],...
        'mod', []...
        );    
    Photometry.data = repmat(data, s.nChannels, 1); % length = # channels

    %% initialize
    originalSamples = max(cellfun(@length, SessionData.demod(:,1))); % just check first channel
    newSamples = ceil(originalSamples/s.downsample);
    allData = NaN(SessionData.nTrials, newSamples);   
    modData = NaN(SessionData.nTrials, newSamples);  % raw data that has not been demodulated         
   
    for fCh=1:s.nChannels
        if s.constantLength
            % baselineEndTime is relative to start of photometry recording
            baselineEndTime = round(SessionData.RawEvents.Trial{1}.States.(s.zeroField)(1,1) - SessionData.RawEvents.Trial{1}.States.(s.startField)(1,1));
        else
            disp('*** variable trial length not yet implemented ***');
        end
        for trial = 1:SessionData.nTrials