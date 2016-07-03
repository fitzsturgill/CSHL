function [avg avgX avgSEM avgAll] = phSessionAvg(SessionData, type, outcome, ax, linespec)

    % I should decimate everything right from the get go at data import???
% 
    fCh = 1; % fluor channel
    sample_rate = 6100;
    dmf = 10;  % decimation factor (must be a factor of sample_rate)
    
    sample_rate = sample_rate / dmf;
    
    trials = bpFilterTrials(SessionData, type, outcome);
%     zeroField = 'DeliverStimulus';
    zeroField = 'PreTrialRecording';
%     startField = 'ITI';
%     endField = 'PostTrialRecording';
    dataLim = [-2, 6];
    % initialize NaN matrix and xdata
    sessionData = NaN(length(trials), round((dataLim(2) - dataLim(1)) * sample_rate));
    avgX = linspace(dataLim(1), dataLim(2), size(sessionData, 2));
    
    %gather trials and align to zero
    for counter = 1:length(trials)
        trial = trials(counter);
        zeroTime = SessionData.RawEvents.Trial{trial}.States.(zeroField)(1,1); % always use first instance and start of zeroField for zeroing (each state has a start and end timestamp
        trialData = SessionData.NidaqData{trial}(:, fCh)';
        trialData = decimate(trialData, dmf); % decimate
%         trialData = smooth(trialData, 10, 'lowess')'; % smooth
        i1 = round((zeroTime + dataLim(1)) * sample_rate);
        % I'm only dealing with trials that don't have enough data after
        % zerotime, not ones that begin too close to zero time (for now):
        i2 = min(i1 + size(sessionData, 2) - 1, size(trialData, 2));        
        sessionData(counter, 1:(i2-i1+1)) = trialData(1, i1:i2);
    end
    
    % convert to deltaF/F
    baselines = nanmean(sessionData(:, 1: abs(dataLim(1)) * sample_rate), 2);
    baselines = repmat(baselines, 1, size(sessionData, 2));
    sessionData = sessionData - baselines;
    sessionData = sessionData ./ baselines;
    
    
    % avg
    avg = nanmean(sessionData);
    avgSEM = nanSEM(sessionData);
    % kludge
    avgAll = sessionData;

    % only plot data if given axis
    if nargin >  3   
        boundedline(avgX, avg, avgSEM, linespec, ax, 'alpha');
    end
        
        
        
        
        