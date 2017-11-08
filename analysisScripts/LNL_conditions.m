% LNL_Conditions

% script to generate trial lookups by condition for
% LNL_odor_v2

%% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks    blocks 2 and 3
    validTrials = filterTE(TE, 'reject', 0);
    Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
    Odor2Trials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);    
    
    rewardTrials = filterTE(TE, 'trialType', 1, 'reject', 0);
    hitTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    missTrials = filterTE(TE, 'trialOutcome', -1, 'reject', 0);
    punishTrials = filterTE(TE, 'trialType', 3, 'reject', 0);    
    neutralTrials = filterTE(TE, 'trialType', [2 4], 'reject', 0);
    block2Trials = filterTE(TE, 'BlockNumber', 2, 'reject', 0);
    block3Trials = filterTE(TE, 'BlockNumber', 3, 'reject', 0);
    csPlusTrials = filterTE(TE, 'trialType', [1 2], 'reject', 0);
    csMinusTrials = filterTE(TE, 'trialType', [3 4], 'reject', 0);    
    trialTypes = 1:4;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end

trialCount = [1:length(TE.filename)]';