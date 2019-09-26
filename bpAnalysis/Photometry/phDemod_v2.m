function demod = phDemod_v2(rawData, refData, refChannel, sampleRate, varargin)



%% optional parameters, first set defaults      
    defaults = {...
        'lowCutoff', 20;... % 15 before
        'forceAmp', 0;... % force demodulation
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    % find channel index
    chix = find(refData.channelsOn == refChannel);
    

    nSamples = length(rawData);
    dt = 1/sampleRate;    
    t = (0:dt:(nSamples - 1) * dt);
    t = t(:);
    

    refData_0 = sin(2*pi*refData.freq(chix)*t + refData.phaseShift(chix));
    refData_90 = sin(2*pi*refData.freq(chix)*t + refData.phaseShift(chix) + pi/2);
    processedData_0 = rawData .* refData_0;
    processedData_90 = rawData .* refData_90;
    %% try filtering first
    % note-   5 pole Butterworth filter in Matlab used in Frohlich and McCormick  
     % Create butterworth filter
    lowCutoff = s.lowCutoff/(sampleRate/2); % filter cutoff normalized to nyquist frequency     
    [z,p,k] = butter(10, lowCutoff, 'low');
    [sos, g] = zp2sos(z,p,k);
    pad = 1; % kludge, always pad
    if pad
        if refData.freq(chix)
            paddedData_0 = processedData_0(1:min(size(processedData_0, 1),sampleRate), 1); % AGV sez: pad with 1s of data, should be in phase as period should be an integer factor of 1 second
            paddedData_90 = processedData_90(1:min(size(processedData_90, 1),sampleRate), 1); % AGV sez: pad with 1s of data, should be in phase as period should be an integer factor of 1 second        
            % HOWEVER- an additional problem is that there is a hardware onset
            %% 
            % transient when the LED turns on
            demodDataFilt_0 = filtfilt(sos,g,[paddedData_0; processedData_0]);
            demodDataFilt_90 = filtfilt(sos,g,[paddedData_90; processedData_90]);                     
            demod_0 = demodDataFilt_0(length(paddedData_0) + 1:end, 1);
            demod_90 = demodDataFilt_90(length(paddedData_90)+1:end, 1);        
            demod = (demod_0 .^2 + demod_90.^2) .^(1/2); % quadrature decoding
        else
            paddedData = rawData(1:min(sampleRate, size(rawData, 1)), 1);
            demod = filtfilt(sos, g, [paddedData; rawData]);
            demod = demod(length(paddedData) + 1:end, 1);
            return
        end            
    end   
    
    % correct for amplitude of DEMODULATION reference     
    % Vsig = Vsig*Vref/2 + Vsig*Vref/2 * Cos(2*Fmod * time)
    % you filter out the second term
    % multiply by two and divide by Vref to get Vsig
    % demod = demod * 2 / amp;
    % amp = 1 for reference so it doesn't matter...

    demod = demod * 2;
    
%     % HOWEVER you still need to get rid of amplitude of MODULATION reference (e.g.
%     % if you multiply a DC offset of 3 by a sinusoid of amplitude 5, you
%     % get a sinusoid of ampltitude 15, so you have to divide by 5 to
%     % recover the signal.
% %     measure amplitude of reference    
%     L = length(rawData);
%     n = 2 ^ nextpow2(L);
%     Y = fft(rawData, n);
%     P2 = abs(Y/n);
%     P1 = P2(1:n/2 + 1);
%     P1(2:end-1) = 2 * P1(2:end-1);
%     f = sampleRate * (0:(n/2))/n;
%     % find indices of mod frequency
%     fix = find(f >= refData.freq(chix), 1);
%     fix = [fix - 1 fix fix + 1]; 
%     amp = max(P1(fix));
% 
% 
%     demod = demod * 2 / amp;