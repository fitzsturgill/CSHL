
DB = dbLoadExperiment('FrankenLNL_RewardPunish');


for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    nTrials = length(TE.filename);
    if ~success
        disp('wtf');
        continue
    end    
    channels = TE.Photometry.settings.channels;
    
    % make sure that Shock times are added to TE.Us (some earlier TEs may not
    % include shock time in TE.Us
    mycount=0;
    if isfield(TE, 'Shock') 
        for trial = 1:nTrials
            if isfinite(TE.Shock{trial}(1))
                TE.Us(trial) = TE.Shock(trial);
            end
        end
    end
    
            
    csWindow = zeros(nTrials, 2);
    csWindow(:,2) = cellfun(@(x,y) y(end) - x(1), TE.Cue2, TE.Trace2);
    csWindowPhasic = [0.3 1.3];

    usWindow = [0 1];  
    
    % percentile value for peak estimations
    percentValue = 0.9;
    
    % estimate respones for different events in each trial for photometry
    % (bpCalcPeak_dFF) and for licking (countEventFromTE)
    TE.licks_cs = countEventFromTE(TE, 'Port1In', csWindow, TE.Cue2);
    TE.lickIntervals_cs = eventIntervalsFromTE(TE, 'Port1In', csWindow, TE.Cue2);
    TE.lickLatency_cs = calcEventLatency(TE, 'Port1In', TE.Cue2, TE.Outcome);
    TE.licks_us = countEventFromTE(TE, 'Port1In', usWindow, TE.Us);
    TE.licks_baseline = countEventFromTE(TE, 'Port1In', [0 4], TE.PreCsRecording);
    TE.lickIntervals_us = eventIntervalsFromTE(TE, 'Port1In', usWindow, TE.Outcome);
    TE.lickLatency_us = calcEventLatency(TE, 'Port1In', TE.Outcome);
    
    for channel = channels
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue2, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue2, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');          
    end     
    
    

    
    
    dbSaveAnimal(DB, animal);            
end