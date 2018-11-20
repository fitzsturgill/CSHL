% SO_RewardPunish_odor_conditions    validTrials = filterTE(TE, 'reject', 0);
    rewardOdorTrials = filterTE(TE, 'trialType', 1:2, 'reject', 0);
    punishOdorTrials = filterTE(TE, 'trialType', 3:4, 'reject', 0);
    uncuedTrials = filterTE(TE, 'trialType', 5:7, 'reject', 0);    
    rewardTrials = filterTE(TE, 'trialOutcome', [1 5], 'reject', 0);
    punishTrials = filterTE(TE, 'trialOutcome', [3 6], 'reject', 0);    
    omitTrials = filterTE(TE, 'trialOutcome', [2 4 7], 'reject', 0);
    trialTypes = 1:7;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end
