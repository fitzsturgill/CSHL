function demod = phDemod_v2(rawData, refData, refChannel, sampleRate, varargin)



%% optional parameters, first set defaults      
    defaults = {...
        'lowCutoff', 15;...
        'forceAmp', 0;... % force demodulation
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    % find channel index
    chix = find(refData.channelsOn == refChannel);
    

    nSamples = length(rawData);
    dt = 1/sampleRate;    
    t = (0:dt:(nSamples - 1) * dt);
    t = t(:);
    
    if s.forceAmp && refData.amp(chix) == 0
        amp = s.forceAmp; % demodulate using nominal LED amplitude even if reference LED is off
    else
        amp = refData.amp(chix);
    end
    refData_0 = sin(2*pi*refData.freq(chix)*t + refData.phaseShift(chix))  * amp;
    refData_90 = sin(2*pi*refData.freq(chix)*t + refData.phaseShift(chix) + pi/2)  * amp;
    processedData_0 = rawData .* refData_0;
    processedData_90 = rawData .* refData_90;
    %% try filtering first
    % note-   5 pole Butterworth filter in Matlab used in Frohlich and McCormick  
     % Create butterworth filter
    lowCutoff = s.lowCutoff/sampleRate * 2; % multiply by 2 to convert to rad/sample- see butter documentation
    % for a cutoff freq of 300Hz and sample rate of 1000Hz, cutoff
    % corresponds to 0.6pi rad/sample    300/1000 * 2 = 0.6    
%     [b, a] = butter(5, lowCutoff, 'low');   
    [z,p,k] = butter(5, lowCutoff, 'low');
    [sos, g] = zp2sos(z,p,k);
    pad = 1;
    if pad
        paddedData_0 = processedData_0(1:sampleRate, 1); % AGV sez: pad with 1s of data, should be in phase as period should be an integer factor of 1 second
        paddedData_90 = processedData_90(1:sampleRate, 1); % AGV sez: pad with 1s of data, should be in phase as period should be an integer factor of 1 second        
        % HOWEVER- an additional problem is that there is a hardware onset
        % transient when the LED turns on
        demodDataFilt_0 = filtfilt(sos,g,[paddedData_0; processedData_0]);
        demodDataFilt_90 = filtfilt(sos,g,[paddedData_90; processedData_90]);        
%         demodDataFilt_0 = filtfilt(b,a,[paddedData_0; processedData_0]);
%         demodDataFilt_90 = filtfilt(b,a,[paddedData_90; processedData_90]);                
        demod_0 = demodDataFilt_0(length(paddedData_0) + 1:end, 1);
        demod_90 = demodDataFilt_90(length(paddedData_90)+1:end, 1);        
    end   
    demod = (demod_0 .^2 + demod_90.^2) .^(1/2); % quadrature decoding
    % correct for amplitude of reference     
    % Vsig = Vsig*Vref/2 + Vsig*Vref/2 * Cos(2*Fmod * time)
    % you filter out the second term
    % multiply by two and divide by Vref to get Vsig
    demod = demod * 2 / amp;