% LNL_Conditions

% script to generate trial lookups by condition for
% LNL_odor_v2

%% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks    blocks 2 and 3
 % Outcomes ->  -1: miss (CS+ specific), 0: false alarm (Cs- specific), 1: hit (Cs+ specific), 2: correct rejection (Cs- specific) (see TrialTypeOutcomePlot) 
    validTrials = filterTE(TE, 'reject', 0);
    Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
    Odor2Trials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);    
    Odor3Trials = filterTE(TE, 'OdorValveIndex', 3, 'reject', 0);  
    uncuedTrials = filterTE(TE, 'OdorValveIndex', 0, 'reject', 0);
    laserTrials = filterTE(TE, 'TriggerPulsePal', 1, 'reject', 0);
    rewardTrials = filterTE(TE, 'ReinforcementOutcome', 'Reward', 'reject', 0);
    hitTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    missTrials = filterTE(TE, 'trialOutcome', -1, 'reject', 0);
    FATrials = filterTE(TE, 'trialOutcome', 0, 'reject', 0); % "false alarm" trials (even if pavlovian)
    CRTrials = filterTE(TE, 'trialOutcome', 2, 'reject', 0); % "correct reject" trials (even if pavlovian)
    punishTrials = filterTE(TE, 'ReinforcementOutcome', 'Punish', 'reject', 0);    
    neutralTrials = filterTE(TE, 'ReinforcementOutcome', 'Neutral', 'reject', 0);
    block2Trials = filterTE(TE, 'BlockNumber', 2, 'reject', 0);
    block3Trials = filterTE(TE, 'BlockNumber', 3, 'reject', 0);
    csPlusTrials = filterTE(TE, 'CSValence', 1, 'reject', 0);
    csMinusTrials = filterTE(TE, 'CSValence', -1, 'reject', 0); 
    uncuedReward = filterTE(TE, 'ReinforcementOutcome', 'Reward', 'OdorValveIndex', 0, 'reject', 0);
    uncuedPunish = filterTE(TE, 'ReinforcementOutcome', 'Punish', 'OdorValveIndex', 0,'reject', 0);
    trialTypes = 1:max(TE.trialType);
    trialsByType = cell(size(trialTypes));
    for ttc = 1:max(TE.trialType)
        trialsByType{ttc} = filterTE(TE, 'trialType', trialTypes(ttc), 'reject', 0);
    end
    trialCount = [1:length(TE.filename)]';
    clear ttc;