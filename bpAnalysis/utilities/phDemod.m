function demod = phDemod(rawData, refData, sampleRate, modRate, lowCutoff)
% DEMODULATE AM-MODULATED INPUT IN QUADRATURE GIVEN A REFERENCE
% tBaseline-  baseline period, in seconds
% lowCutoff corner frequency for 5-pole butterworth filter (lowpass),
% default = [], i.e. no filtering is performed

    if nargin < 5
        lowCutoff = []; 
    end
        
%     if nargin < 5
%         tBaseline = []; % if empty, normalize (zscore) by entire range
%     end
    
    if size(rawData, 2) ~= 1 ||size(refData, 2) ~= 1
        disp('*** Error in phDemod, refData and rawData must be column vectors ***');
        demod = [];
        return
    end
%     rawData(:,1) = rawData; % ensure column vectors
%     refData(:,1) = refData; 


    nSamples = length(rawData);

    refData = refData(1:nSamples,1); % shorten refData to same size as rawData
    
    refData = refData - mean(refData); % *** get rid of DC offset!!!!
    
%     if ~isempty(tBaseline)
%         useForBaseline = rawData(1:floor(tBaseline / (1/sampleRate)));        
% %         [processedData, ~, ~] = zscoreByRange(rawData, 1, preSamples);
%     else
%         useForBaseline = rawData; 
% %         [processedData, ~, ~] = zscore(rawData);        % just use the whole thing
%     end
%     blMean = mean(useForBaseline);
%     blSD = std(useForBaseline);    
%     refData = 
%     [refData, ~, ~] = zscore(refData);

    % generate 90degree shifted copy of refData
    samplesPerPeriod = 1/modRate / (1/sampleRate);
    quarterPeriod = round(samplesPerPeriod / 4); % ideally you shouldn't have to round, i.e. mod frequencies should be close to factors of sample freq
    refData90 = circshift(refData, [1 quarterPeriod]);

    processedData_0 = rawData .* refData;
    processedData_90 = rawData .* refData90;

    demodData = (processedData_0 .^2 + processedData_90 .^2) .^(1/2); % quadrature decoding

     % Create butterworth filter
    lowCutoff = lowCutoff/sampleRate * 2; % multiply by 2 to convert to rad/sample- see butter documentation
    % for a cutoff freq of 300Hz and sample rate of 1000Hz, cutoff
    % corresponds to 0.6pi rad/sample    300/1000 * 2 = 0.6

    % note-   5 pole Butterworth filter in Matlab used in Frohlich and McCormick  
    [b, a] = butter(5, lowCutoff, 'low');   
    pad = 1;
    if pad
%         paddedData = fliplr(demodData(1:sampleRate, 1)); % pad with 1s of reflected data
%         paddedData = demodData(randperm(sampleRate), 1); % pad with 1s of randomized data (should still contain DC trend)
        paddedData = demodData(1:sampleRate, 1); % AGV sez: pad with 1s of data, should be in phase as period should be an integer factor of 1 second
        % HOWEVER- an additional problem is that there is a hardware onset
        % transient when the LED turns on
        demodDataFilt = filtfilt(b,a,[paddedData; demodData]);
        demod = demodDataFilt(length(paddedData) + 1: end, 1);
    else
        demod = filtfilt(b, a, demodData);
    end
    
    % correct for amplitude of reference 
    
    % Vsig = Vsig*Vref/2 + Vsig*Vref/2 * Cos(2*Fmod * time)
    % you filter out the second term
    % multiply by two and divide by Vref to get Vsig
    
    modAmp = calcSinusoidAmp(refData);
    demod = demod * 2 / modAmp;
%     fig = ensureFigure('test', 1);
%     plot(demodDataFilt);

