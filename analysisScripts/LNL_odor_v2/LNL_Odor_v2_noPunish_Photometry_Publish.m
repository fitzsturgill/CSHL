
% function TE = LNL_Odor_v2_noPunish_Photometry_Publish(DB, animal, reloadPhotometry, 
% LNL_Odor_v2_noPunish_Photometry_Publish

DB = dbLoadExperiment('reversals_noPunish_publish');

% deconvolve photometry data
tau = 2; % time constant of exponential decay kernel
kDuration = 3; % duration of kernel
Fs = 20;
kt = (1:kDuration*20)*(1/Fs) - (1/Fs); 
k = exp(-1 * (1/tau) * kt);
k = k / trapz(k);
% deconvolve data
epsilon = 0.1;
PhotometryFields = {'Photometry', 'PhotometryExpFit'};
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    nTrials = length(TE.filename);
    if ~success
        disp('wtf');
        continue
    end    

% deconvolve photometry data
    for pCounter = 1:length(PhotometryFields)
        pField = PhotometryFields{pCounter};
        deconvSettings = struct('epsilon', [], 'tau', [], 'kernelLength', []);
        for channel = 1:2
            TE.(pField).data(channel).dFdeconv = bpDeconv(TE.(pField).data(channel).dF, k, epsilon, 'none');
            deconvSettings.epsilon(end+1) = epsilon;
            deconvSettings.tau(end+1) = tau;
            deconvSettings.kernelLength(end+1) = kDuration;
            %         recompute the deconvolved z-score
            TE.(pField).data(channel).ZSdeconv = TE.(pField).data(channel).dFdeconv ./ nanmean(nanstd(TE.(pField).data(channel).dFdeconv(:,bpX2pnt(1, 20):bpX2pnt(4, 20)), 0, 2));
        end
        TE.(pField).settings.deconvSettings = deconvSettings;
    end
    
    csWindow = zeros(nTrials, 2);
    csWindow(:,2) = cellfun(@(x,y,z) max(x(end), y(end)) - z(1), TE.AnswerLick, TE.AnswerNoLick, TE.Cue);
    csWindowPhasic = [0.3 1.3];
    
    usWindow = [0 1];   
    
    % percentile value for peak estimations
    percentValue = 0.8;
    
    % estimate respones for different events in each trial for photometry
    % (bpCalcPeak_dFF) and for licking (countEventFromTE)
    TE.licks_cs = countEventFromTE(TE, 'Port1In', csWindow, TE.Cue);
    TE.licks_us = countEventFromTE(TE, 'Port1In', [0 2], TE.Us);
    TE.licks_baseline = countEventFromTE(TE, 'Port1In', [0 4], TE.PreCsRecording);
    
    pupLag = 0.3;
    TE.pupil_cs = bpCalcPeak_Pupil(TE.pupil, 'zeroTimes', TE.Cue, 'window', csWindow);
    TE.pupil_us = bpCalcPeak_Pupil(TE.pupil, 'zeroTimes', TE.Us, 'window', usWindow + pupLag);
    TE.pupil_baseline = bpCalcPeak_Pupil(TE.pupil, 'zeroTimes', TE.PreCsRecording, 'window', [0 4]);
        
    TE.whisk_cs = bpCalcPeak_Whisk(TE.Whisk, 'zeroTimes', TE.Cue, 'window', csWindow);
    TE.whisk_us = bpCalcPeak_Whisk(TE.Whisk, 'zeroTimes', TE.Us, 'window', usWindow);
    TE.whisk_baseline = bpCalcPeak_Whisk(TE.Whisk, 'zeroTimes', TE.PreCsRecording, 'window', [0 4]);    
                
    TE.wheel_cs = bpCalcPeak_Wheel(TE.Wheel, 'zeroTimes', TE.Cue, 'window', csWindow);
    TE.wheel_us = bpCalcPeak_Wheel(TE.Wheel, 'zeroTimes', TE.Us, 'window', usWindow);
    TE.wheel_baseline = bpCalcPeak_Wheel(TE.Wheel, 'zeroTimes', TE.PreCsRecording, 'window', [0 4]);    
    
    for channel = channels                
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
    end 
    
    dbSaveAnimal(DB, animal);    
end