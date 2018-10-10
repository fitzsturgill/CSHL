% reversals_noPunish_poolReversals

DB = dbLoadExperiment('reversals_noPunish_publish');

for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    if ~success
        disp('wtf');
        continue
    end    
    
    
    dataToPull = {...
        'licks_cs', TE.licks_cs.rate,...
        'licks_us', TE.licks_us.rate,...
        'licks_baseline', TE.licks_baseline.rate,...
        'pupil_cs', TE.pupil_cs.data,...
        'pupil_csBaselined', TE.pupil_csBaselined.data,...        
        'pupil_us', TE.pupil_us.data,...
        'pupil_baseline', TE.pupil_baseline.data,... 
        'whisk_cs', TE.whisk_cs.data,...
        'whisk_us', TE.whisk_us.data,...
        'whisk_baseline', TE.whisk_baseline.data,...                
        'wheel_cs', TE.wheel_cs.data,...
        'wheel_us', TE.wheel_us.data,...
        'wheel_baseline', TE.wheel_baseline.data,...
        'phPeakMean_cs_ch1', TE.phPeakMean_cs(1).data,...
        'phPeakPercentile_cs_ch1', TE.phPeakPercentile_cs(1).data,...        
        'phPeakMean_us_ch1', TE.phPeakMean_us(1).data,...
        'phPeakPercentile_us_ch1', TE.phPeakPercentile_us(1).data,...        
        'phBaseline_ch1', TE.phPeakMean_baseline(1).data,...
        'phPeakMean_cs_ch2', TE.phPeakMean_cs(2).data,...
        'phPeakPercentile_cs_ch2', TE.phPeakPercentile_cs(2).data,...        
        'phPeakMean_us_ch2', TE.phPeakMean_us(2).data,...
        'phPeakPercentile_us_ch2', TE.phPeakPercentile_us(2).data,...        
        'phBaseline_ch2', TE.phPeakMean_baseline(2).data,...        
        'trialOutcome', TE.trialOutcome,...
        'trialType', TE.trialType,...
        'trialNumber', TE.trialNumber,...
        'filename', TE.filename,...
        'ReinforcementOutcome', TE.ReinforcementOutcome,...
        'OdorValveIndex', TE.OdorValveIndex,...
        'csLicksROC', TE.AnswerLicksROC,...
        };
    
    RE = struct();
    RE.csPlus = extractReversalsFromTE(TE, csPlusTrials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    RE.csMinus = extractReversalsFromTE(TE, csMinusTrials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    RE.thirdOdor = extractReversalsFromTE(TE, Odor3Trials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    RE.csPlusReward = extractReversalsFromTE(TE, csPlusTrials & rewardTrials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    RE.allTrials = extractReversalsFromTE(TE, validTrials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    

    save(fullfile(DB.path, 'pooled', ['RE_' animal '.mat']), 'RE');
    disp(['*** saving: ' fullfile(DB.path, 'pooled', ['RE_' animal '.mat']) ' ***']);
end
    
