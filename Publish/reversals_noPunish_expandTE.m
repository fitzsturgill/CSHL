
% function TE = LNL_Odor_v2_noPunish_Photometry_Publish(DB, animal, reloadPhotometry, 
% LNL_Odor_v2_noPunish_Photometry_Publish

DB = dbLoadExperiment('reversals_noPunish_publish');
channels = [1 2];
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
% for counter = 1:length(DB.animals)
for counter = 7
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
    TE.lickIntervals_cs = eventIntervalsFromTE(TE, 'Port1In', csWindow, TE.Cue);
    TE.licks_us = countEventFromTE(TE, 'Port1In', [0 2], TE.Us);
    TE.licks_baseline = countEventFromTE(TE, 'Port1In', [0 4], TE.PreCsRecording);
    
    pupLag = 0.3; pupilWindow = [0.5 csWindow(2)]; % pupil never really ascends until nearly a second after the cue turns on
    TE.pupil_cs = bpCalcPeak_Pupil(TE.pupil, 'zeroTimes', TE.Cue, 'window', csWindow + pupLag);
    TE.pupil_us = bpCalcPeak_Pupil(TE.pupil, 'zeroTimes', TE.Us, 'window', usWindow + pupLag);
    TE.pupil_baseline = bpCalcPeak_Pupil(TE.pupil, 'zeroTimes', TE.PreCsRecording, 'window', [0 4]);
    TE.pupil_csBaselined.data = TE.pupil_cs.data - TE.pupil_baseline.data;
        
    TE.whisk_cs = bpCalcPeak_Whisk(TE.Whisk, 'zeroTimes', TE.Cue, 'window', csWindow);
    TE.whisk_us = bpCalcPeak_Whisk(TE.Whisk, 'zeroTimes', TE.Us, 'window', usWindow);
    TE.whisk_baseline = bpCalcPeak_Whisk(TE.Whisk, 'zeroTimes', TE.PreCsRecording, 'window', [0 4]);    
    TE.whisk_csBaselined.data = TE.whisk_cs.data - TE.whisk_baseline.data;                

    TE.wheel_cs = bpCalcPeak_Wheel(TE.Wheel, 'zeroTimes', TE.Cue, 'window', csWindow);
    TE.wheel_us = bpCalcPeak_Wheel(TE.Wheel, 'zeroTimes', TE.Us, 'window', usWindow);
    TE.wheel_baseline = bpCalcPeak_Wheel(TE.Wheel, 'zeroTimes', TE.PreCsRecording, 'window', [0 4]);    
    

    for channel = channels
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
        
        TE.phPeakMean_cs_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'mean', 'phField', 'ZSdeconv');
        TE.phPeakMean_us_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'mean', 'phField', 'ZSdeconv');
        TE.phPeakMean_baseline_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZSdeconv');
        TE.phPeakPercentile_cs_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZSdeconv');
        TE.phPeakPercentile_us_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZSdeconv');
    end 
    
    %% find best window for cue response, by discriminability of CS+ and CS-
    % responses ("hit" vs "correct rejection")
%     spacing = 0.1; % 0.1s steps to determine optimum window
    assert(Fs == 20, 'Sample rate assumed to be 20');
    commonCueWindow = [0 min(csWindow(:,2))]; % in case you've changed the delay across sessions, base windows upon smallest delay
%     find peak us response
    avgData = phAverageFromTE(TE, csPlusTrials & hitTrials, 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue, 'window', commonCueWindow); %high value, reward
%     ensureFigure('test', 1);
%     plot(avgData.xData, avgData.Avg);
    [~, mix] = max(avgData.Avg);
%     mix = fix(mix/2)*2; % make even
    prePoints = (1:2:mix)';
    postPoints = mix+1:2:bpX2pnt(commonCueWindow(2), 20, 0);
    % make matrix of pre, post points, and D' for each pairing
    preMatrix = zeros(length(prePoints), length(postPoints));
    postMatrix = zeros(length(prePoints), length(postPoints));
    dMatrix = zeros(length(prePoints), length(postPoints));
    
    % collect phWindow for commonCuewindow for each channel
    for ch = 1:2        
        [phData_plus, ~] = phAlignedWindow(TE, csPlusTrials & hitTrials, ch, 'zeroTimes', TE.Cue, 'window', commonCueWindow);
        [phData_minus, ~] = phAlignedWindow(TE, csMinusTrials & CRTrials, ch, 'zeroTimes', TE.Cue, 'window', commonCueWindow);
        for counter = 1:numel(dMatrix)
            [D, P] = rocarea(nanmean(phData_plus(:,preMatrix(counter):postMatrix(counter)), 2), nanmean(phData_plus(:,preMatrix(counter):postMatrix(counter)), 2), 'scale');
            dMatrix(counter) = D;
        end
    end
        
        
    

    
    
    
    
    
%%     dbSaveAnimal(DB, animal);        
end