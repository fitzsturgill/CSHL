function demod = phDemod(rawData, refData, sampleRate, modRate, lowCutoff)
% DEMODULATE AM-MODULATED INPUT IN QUADRATURE GIVEN A REFERENCE
% tBaseline-  baseline period, in seconds
% lowCutoff corner frequency for 5-pole butterworth filter (lowpass),
% default = [], i.e. no filtering is performed

    if nargin < 5
        lowCutoff = []; 
    end
    
    if isnumeric(refData)
        demodMode = 1;
        refData = refData(:); % ensure column vector
        rawData = rawData(:); % ensure column vector
    elseif isstruct(refData)
        demodMode = 2;
    end

    nSamples = length(rawData);
    switch demodMode
        case 1
            refData_0 = refData(1:nSamples,1); % shorten refData to same size as rawData    
            refData_0 = refData_0 - mean(refData_0); % *** get rid of DC offset!!!! 
            % generate 90degree shifted copy of refData
            samplesPerPeriod = 1/modRate / (1/sampleRate);
            quarterPeriod = round(samplesPerPeriod / 4); % ideally you shouldn't have to round, i.e. mod frequencies should be close to factors of sample freq
            refData_90 = circshift(refData_0, [1 quarterPeriod]);
        case 2
            dt = 1/sampleRate;    
            t = (0:dt:(nSamples - 1) * dt);
            t = t(:);
            refData = sin(2*pi*refData.freq*t + refData.)  * S.GUI.LED1_amp;

    processedData_0 = rawData .* refData_0;
    processedData_90 = rawData .* refData_90;
    %% try filtering first
    % note-   5 pole Butterworth filter in Matlab used in Frohlich and McCormick  
     % Create butterworth filter
    lowCutoff = lowCutoff/sampleRate * 2; % multiply by 2 to convert to rad/sample- see butter documentation
    % for a cutoff freq of 300Hz and sample rate of 1000Hz, cutoff
    % corresponds to 0.6pi rad/sample    300/1000 * 2 = 0.6    
    [b, a] = butter(5, lowCutoff, 'low');   
    pad = 1;
    if pad
%         paddedData = fliplr(demodData(1:sampleRate, 1)); % pad with 1s of reflected data
%         paddedData = demodData(randperm(sampleRate), 1); % pad with 1s of randomized data (should still contain DC trend)
        paddedData_0 = processedData_0(1:sampleRate, 1); % AGV sez: pad with 1s of data, should be in phase as period should be an integer factor of 1 second
        paddedData_90 = processedData_90(1:sampleRate, 1); % AGV sez: pad with 1s of data, should be in phase as period should be an integer factor of 1 second        
        % HOWEVER- an additional problem is that there is a hardware onset
        % transient when the LED turns on
        demodDataFilt_0 = filtfilt(b,a,[paddedData_0; processedData_0]);
        demodDataFilt_90 = filtfilt(b,a,[paddedData_90; processedData_90]);        
        demod_0 = demodDataFilt_0(length(paddedData_0) + 1:end, 1);
        demod_90 = demodDataFilt_90(length(paddedData_90)+1:end, 1);        
    else
%         demod_0 = filtfilt(b, a, demodData_0);
    end
    
    
    
    demod = (demod_0 .^2 + demod_90.^2) .^(1/2); % quadrature decoding




    
    % correct for amplitude of reference 
    
    % Vsig = Vsig*Vref/2 + Vsig*Vref/2 * Cos(2*Fmod * time)
    % you filter out the second term
    % multiply by two and divide by Vref to get Vsig

    modAmp = calcSinusoidAmp(refData_0);
    demod = demod * 2 / modAmp;
%     fig = ensureFigure('test', 1);
%     plot(demodDataFilt);