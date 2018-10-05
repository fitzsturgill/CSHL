
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
    
    usWindow = [0 0.75];
    usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'
    
    % percentile value for peak estimations
    percentValue = 0.8;
    
    % estimate respones for different events in each trial for photometry
    % (bpCalcPeak_dFF) and for licking (countEventFromTE)
    TE.csLicks = countEventFromTE(TE, 'Port1In', [0 3], TE.Cue);
    TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);

    for channel = channels
        TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
        TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'mean', 'phField', 'ZS');
        if channel == 1
            TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
            TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
            TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
        elseif channel == 2
            TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
            TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');        
            TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');  
        end
    end    
end