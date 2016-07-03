function [avg, avgX, avgSEM] = demodSession_old(SessionData, type, outcome, ax, linespec)

    % I should decimate everything right from the get go at data import???
% 
    fCh = 1; % fluor channel
    sampleRate = 6100;
    modF = 211;    
    dmf = 1;  % decimation factor (must be a factor of sample_rate)
    
    sampleRate = sampleRate / dmf;
    
    trials = bpFilterTrials(SessionData, type, outcome);
    zeroField = 'DeliverStimulus';
%     startField = 'ITI';
%     endField = 'PostTrialRecording';
    dataLim = [-2, 4];
    % initialize NaN matrix and xdata

    avgX = linspace(dataLim(1), dataLim(2), size(SessionData, 2));
    
    %gather trials and align to zero
    for counter = 1:length(trials)
        trial = trials(counter);
        zeroTime = SessionData.RawEvents.Trial{trial}.States.(zeroField)(1,1); % always use first instance and start of zeroField for zeroing (each state has a start and end timestamp
        rawData = SessionData.NidaqData{trial,1}(:,1);
        nSamples = size(rawData, 1);
        refData = SessionData.NidaqData{trial,2}(1:nSamples,1);
        [finalData, dmf, dmp] = decode_lockin_fn(rawData, refData, [], modF, sampleRate, 0);
        
%         trialData = SessionData.NidaqData{trial}(:, fCh)';
%         trialData = decimate(trialData, dmf); % decimate
%         trialData = smooth(trialData, 50, 'lowess')'; % smooth
%         i1 = round((zeroTime + dataLim(1)) * sample_rate);
        % I'm only dealing with trials that don't have enough data after
        % zerotime, not ones that begin too close to zero time (for now):
%         i2 = min(i1 + size(SessionData, 2) - 1, size(trialData, 2));   
        if counter == 1
            allData = NaN(length(trials), size(finalData, 1));
            avgX = linspace(0, nSamples/sampleRate, nSamples);
        end
        allData(counter, :) = finalData;
    end
    
    % convert to deltaF/F
%     baselines = nanmean(SessionData(:, 1: abs(dataLim(1)) * sample_rate), 2);
%     baselines = repmat(baselines, 1, size(SessionData, 2));
%     SessionData = SessionData - baselines;
%     SessionData = SessionData ./ baselines;
%     
    
    % avg
    avg = nanmean(allData);
    avgSEM = std(allData, 'omitnan') ./ sqrt(sum(~isnan(allData), 1));

    % only plot data if given axis
    if nargin >  3   
        boundedline(avgX, avg, avgSEM, linespec, ax, 'alpha');
    end
        
        
        
        
        