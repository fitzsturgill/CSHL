% SO_varyRewardSize_odorV2_conditions    validTrials = filterTE(TE, 'reject', 0);
    cuedTrials = filterTE(TE, 'trialType', 1:3, 'reject', 0);   
    uncuedTrials = filterTE(TE, 'trialType', 4:6, 'reject', 0);    
    smallRewardTrials = filterTE(TE, 'rewardSize', 2, 'reject', 0);
    bigRewardTrials = filterTE(TE, 'rewardSize', 8, 'reject', 0);    
    omitTrials = filterTE(TE, 'trialType', [3 6], 'reject', 0);
    trialTypes = 1:6;
    trialsByType = cell(size(trialTypes));
    for sss = 1:length(trialTypes)
        trialsByType{sss} = filterTE(TE, 'trialType', trialTypes(sss), 'reject', 0);
    end
    
    clear sss