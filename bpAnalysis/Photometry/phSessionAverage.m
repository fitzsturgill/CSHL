function [avg, avgX, avgSEM] = phSessionAverage(SessionData, type, outcome, fCh, zeroField)
    %% assumes consistent timing of events relative to photometry acquisition

    if nargin < 5
        zeroField = 'DeliverStimulus';
    end
    trials = bpFilterTrials(SessionData, type, outcome);
    
    try
        sampleRate = SessionData.Settings.sample_rate;
    catch
        sampleRate = 6100;
    end    


% initialize
    allData = NaN(length(trials), size(SessionData.demod{trials(1), fCh}, 1));    
    zeroTime = round(SessionData.RawEvents.Trial{trials(1)}.States.(zeroField)(1,1) - SessionData.RawEvents.Trial{trials(1)}.States.PreTrialRecording(1,1));
    % why isn't the above line consistently working????
%     zeroTime = 2;
    blSamples = zeroTime * sampleRate;
    
    avgX = linspace(-zeroTime, (size(allData, 2) - 1) / sampleRate, size(allData, 2));

    
    for counter = 1:length(trials)
        trial = trials(counter);
        trialData = SessionData.demod{trial, fCh};        
        trialData = trialData'; % convert to row vector
        allData(counter, :) = trialData(1, 1:size(allData, 2));
    end
    
    % convert to deltaF/F    
    % avg, use of nanmean not currently necessary, may be so in future
    blF = nanmean(allData(:, 1:blSamples));
    allData = allData / mean(blF);
    avg = nanmean(allData);
    avgSEM = std(allData, 'omitnan') ./ sqrt(sum(~isnan(allData), 1));

        
        
        