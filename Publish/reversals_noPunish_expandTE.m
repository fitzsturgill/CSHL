
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
for counter = 1:7
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

    saveName = sprintf('LNL_reversals_optimize_windows_Cue_%s', animal);     
    fhc(end + 1) = ensureFigure(saveName, 1);
    mcLandscapeFigSetup(gcf);
    
    assert(Fs == 20, 'Sample rate assumed to be 20');
    commonCueWindow = [0 min(csWindow(:,2))]; % in case you've changed the delay across sessions, base windows upon smallest delay

    
    % collect phWindow for commonCuewindow for each channel
    
    optPlusTrials = csPlusTrials & hitTrials;
    % substitute bottom quintile of cs licks to account for fact that
    % complete absence of licks (CRTrials) is likely to be an overly
    % sensitive measure of low reward expectation       
    optMinusTrials = csMinusTrials & (TE.licks_cs.count <= percentile(TE.licks_cs.count, 0.2));
    for channel = 1:2  
        
    %     find peak cs response
        avgData = phAverageFromTE(TE, csPlusTrials & hitTrials, channel, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue, 'window', commonCueWindow); %high value, reward
        [mv, mix] = max(avgData.Avg);
        % if there is no peak, use point 2/3 way into window
        if (mix > 0.9 * length(avgData.Avg))
            mix = round(1/3 * length(avgData.Avg));
            mv = avgData.Avg(mix);
        end
        prePoints = (1:2:mix)';
        postPoints = mix+1:2:bpX2pnt(commonCueWindow(2), 20, 0);
        % make matrix of pre, post points, and auROC, scaled for each pairing
        preMatrix = repmat(prePoints, 1, length(postPoints));
        postMatrix = repmat(postPoints, length(prePoints), 1);
        window_auROC = zeros(length(prePoints), length(postPoints));    
        [phData_plus, ~] = phAlignedWindow(TE, optPlusTrials, channel, 'zeroTimes', TE.Cue, 'window', commonCueWindow, 'FluorDataField', 'ZS'); 
        [phData_minus, ~] = phAlignedWindow(TE, optMinusTrials, channel, 'zeroTimes', TE.Cue, 'window', commonCueWindow, 'FluorDataField', 'ZS');
        for pcounter = 1:numel(preMatrix)
            [D, P] = rocarea(nanmean(phData_plus(:,preMatrix(pcounter):postMatrix(pcounter)), 2), nanmean(phData_minus(:,preMatrix(pcounter):postMatrix(pcounter)), 2), 'scale');
            window_auROC(pcounter) = D;
        end
        [m, mwix] = max(window_auROC(:));
        csWindowOptIx = [preMatrix(mwix) postMatrix(mwix)];
        csWindowOpt = [bpPnt2x(csWindowOptIx(1), Fs) bpPnt2x(csWindowOptIx(2), Fs)];
        TE.phPeakMean_optimize_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindowOpt, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_optimize_cs(channel).settings.optimumWindowIx = csWindowOptIx; % in points from Cue onset
        TE.phPeakMean_optimize_cs(channel).settings.optimumWindow = csWindowOpt; % in seconds
        TE.phPeakMean_optimize_cs(channel).settings.optimumWindow_avg = avgData.Avg;
        TE.phPeakMean_optimize_cs(channel).settings.optimumWindow_avgX = avgData.xData;
        TE.phPeakMean_optimize_cs(channel).settings.optimumWindow_watershedPoint = mix;  
        TE.phPeakMean_optimize_cs(channel).settings.optimumWindow_StartPointMatrix = avgData.xData(preMatrix);
        TE.phPeakMean_optimize_cs(channel).settings.optimumWindow_StopPointMatrix = avgData.xData(postMatrix);
        TE.phPeakMean_optimize_cs(channel).settings.optimumWindow_auROC = window_auROC;
        TE.phPeakMean_optimize_cs(channel).settings.optimumWindow_auROCMax = window_auROC(mwix);

        subplot(2,2,channel);  surf(TE.phPeakMean_optimize_cs(channel).settings.optimumWindow_StartPointMatrix,...
            TE.phPeakMean_optimize_cs(channel).settings.optimumWindow_StopPointMatrix,...
            window_auROC); hold on; 
        plot3(avgData.xData(preMatrix(mwix)), avgData.xData(postMatrix(mwix)), window_auROC(mwix), 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
        title(sprintf('%s: channel %u', animal, channel), 'Interpreter', 'none');        
        xlabel('time from cue (s)'); ylabel('time from cue (s)');
        subplot(2,2,channel + 2);
        plot(avgData.xData, avgData.Avg); hold on;
        stem(avgData.xData(mix), avgData.Avg(mix), 'k');
        ylabel('ZS'); xlabel('time from cue (s)');
        addStimulusPatch(gca, [csWindowOpt 0 0.75], sprintf('auROC=%.2f, Nplus=%u, Nminus=%u', window_auROC(mwix),...
            sum(optPlusTrials), sum(optMinusTrials)), [1 0 0]);
        
        if saveOn
            saveas(gcf, fullfile(animalpath, [saveName '.fig']));
            saveas(gcf, fullfile(animalpath, [saveName '.jpg']));   
        end        
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