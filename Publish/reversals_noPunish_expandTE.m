
% function TE = LNL_Odor_v2_noPunish_Photometry_Publish(DB, animal, reloadPhotometry, 
% LNL_Odor_v2_noPunish_Photometry_Publish
saveOn = 1;
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
fhc = [];
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    animalpath = fullfile(DB.path, 'animals', animal, filesep);
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
    
    % smooth lick times into lick rates
    maxLickRate = 12;
    normalLickRate = 8.2;
    window = [-7 4];
    sigma = 0.2; % standard deviation of gaussian to smooth licks
    binEdges = linspace(window(1), window(2), TE.Photometry.sampleRate * diff(window) + 1);
    [eventTimes, eventTrials] = extractEventTimesFromTE(TE, 1:length(TE.filename), 'Port1In', 'zeroTimes', TE.Us, 'window', window);
    counts = histCountsByTrial(eventTimes, eventTrials, binEdges);
    counts = counts .* TE.Photometry.sampleRate;
    kernel_time = -3*sigma:1/TE.Photometry.sampleRate:3*sigma; % kernel is 3 standard deviations wide 
    kernel = normpdf(kernel_time, 0, sigma);
    kernel = kernel / sum(kernel); % normalize kernel
    TE.lickRates.data = conv2(counts, kernel, 'same');
    TE.lickRates.data(TE.lickRates.data > maxLickRate) = normalLickRate;
    TE.lickRates.settings.sigma = sigma;
    TE.lickRates.startTime = cellfun(@(x) x(1), TE.Us) + window(1);
    TE.lickRates.sampleRate = TE.Photometry.sampleRate;
    TE.lickRates.maxLickRate = maxLickRate;
    
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
        
        TE.phPeakMean_cs_expFit(channel) = bpCalcPeak_dFF(TE.PhotometryExpFit, channel, csWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_us_expFit(channel) = bpCalcPeak_dFF(TE.PhotometryExpFit, channel, usWindow, TE.Us, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_baseline_expFit(channel) = bpCalcPeak_dFF(TE.PhotometryExpFit, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs_expFit(channel) = bpCalcPeak_dFF(TE.PhotometryExpFit, channel, csWindow, TE.Cue, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
        TE.phPeakPercentile_us_expFit(channel) = bpCalcPeak_dFF(TE.PhotometryExpFit, channel, usWindow, TE.Us, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
        
        TE.phPeakMean_cs_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'mean', 'phField', 'ZSdeconv');
        TE.phPeakMean_us_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'mean', 'phField', 'ZSdeconv');
        TE.phPeakMean_baseline_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZSdeconv');
        TE.phPeakPercentile_cs_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZSdeconv');
        TE.phPeakPercentile_us_deconv(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZSdeconv');
    end 
    
    %% find best window for cue response, by discriminability of CS+ and CS-
    % responses ("hit" vs "correct rejection")
%     spacing = 0.1; % 0.1s steps to determine optimum window
    minWindow = 0.5;
    maxWindow = 2;
    window = [0 2];
   
    
    
    optPlusTrials = csPlusTrials & hitTrials;
    % substitute bottom quintile of cs licks to account for fact that
    % complete absence of licks (CRTrials) is likely to be an overly
    % sensitive measure of low reward expectation       
    optMinusTrials = csMinusTrials & (TE.licks_cs.count <= percentile(TE.licks_cs.count, 0.2));
    
    % standard Photometry field
    saveName = sprintf('LNL_reversals_optimize_windows_Cue_%s', animal);
    fhc(end + 1) = ensureFigure(saveName, 1);     
    TE.phPeakMean_optimize_cs = bpCalcPeak_optimumWindow_dFF(TE, optPlusTrials, optMinusTrials,...
        'plot', true, 'PhotometryField', 'Photometry', 'window', window, 'zeroTimes', TE.Cue, 'minWindow', minWindow, 'maxWindow', maxWindow, 'fig', gcf);
    title(gca, sprintf('%s: channel %u', animal, channel), 'Interpreter', 'none');        
    if saveOn
        saveas(gcf, fullfile(animalpath, [saveName '.fig']));
        saveas(gcf, fullfile(animalpath, [saveName '.jpg']));   
    end       
    
    % expFit Photometry field
    saveName = sprintf('LNL_reversals_optimize_windows_Cue_expFit_%s', animal);
    fhc(end + 1) = ensureFigure(saveName, 1);     
    TE.phPeakMean_optimize_cs = bpCalcPeak_optimumWindow_dFF(TE, optPlusTrials, optMinusTrials,...
        'plot', true, 'PhotometryField', 'PhotometryExpFit', 'window', window, 'zeroTimes', TE.Cue, 'minWindow', minWindow, 'maxWindow', maxWindow, 'fig', gcf);
    title(gca, sprintf('ExpFit: %s: channel %u', animal, channel), 'Interpreter', 'none');        
    if saveOn
        saveas(gcf, fullfile(animalpath, [saveName '.fig']));
        saveas(gcf, fullfile(animalpath, [saveName '.jpg']));   
    end       
    



    dbSaveAnimal(DB, animal);        
end
h = waitbar(0, 'slowly writing pdfs');

pdfPlus = fullfile(DB.path, 'pooled', 'CueWindow_optmimize_pooled.pdf');
for counter = 1:length(fhc)    
    if counter == 1
        export_fig(fhc(counter),pdfPlus);  % write to pdf
    else
        export_fig(fhc(counter),'-append',pdfPlus);  % write to pdf
    end
    waitbar(counter/length(fhc));
end
close(h);