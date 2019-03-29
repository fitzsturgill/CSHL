% LNL_Conditions

% script to generate trial lookups by condition for
% LNL_odor_v2

%% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks    blocks 2 and 3
 % Outcomes ->  -1: miss (CS+ specific), 0: false alarm (Cs- specific), 1: hit (Cs+ specific), 2: correct rejection (Cs- specific) (see TrialTypeOutcomePlot) 
    validTrials = filterTE(TE, 'reject', 0);
    Odor1Valve1Trials = filterTE(TE, 'Odor1ValveIndex', 1, 'reject', 0);
    Odor1Valve2Trials = filterTE(TE, 'Odor1ValveIndex', 2, 'reject', 0);    
    Odor1Valve3Trials = filterTE(TE, 'Odor1ValveIndex', 3, 'reject', 0);  
    Odor1Valve4Trials = filterTE(TE, 'Odor1ValveIndex', 4, 'reject', 0);      
    Odor1Valve5Trials = filterTE(TE, 'Odor1ValveIndex', 5, 'reject', 0);       
    Odor2Valve1Trials = filterTE(TE, 'Odor2ValveIndex', 1, 'reject', 0);
    Odor2Valve2Trials = filterTE(TE, 'Odor2ValveIndex', 2, 'reject', 0);    
    Odor2Valve3Trials = filterTE(TE, 'Odor2ValveIndex', 3, 'reject', 0);  
    Odor2Valve4Trials = filterTE(TE, 'Odor2ValveIndex', 4, 'reject', 0);   
    Odor2Valve5Trials = filterTE(TE, 'Odor2ValveIndex', 5, 'reject', 0);   
    TinyPuffTrials = filterTE(TE, 'Odor2ValveIndex', -1, 'reject', 0);     
    uncuedTrials = filterTE(TE, 'Odor1ValveIndex', 0, 'Odor2ValveIndex', 0, 'reject', 0);
    rewardTrials = filterTE(TE, 'ReinforcementOutcome', 'Reward', 'reject', 0);
    punishTrials = filterTE(TE, 'ReinforcementOutcome', 'Punish', 'reject', 0);    
    neutralTrials = filterTE(TE, 'ReinforcementOutcome', 'Neutral', 'reject', 0);
    if isfield(TE, 'ShockCurrent')
        shockTrials = filterTE(TE, 'ReinforcementOutcome', 'Shock', 'reject', 0) & (abs(TE.ShockCurrent) > 200);
    else
        shockTrials = filterTE(TE, 'ReinforcementOutcome', 'Shock', 'reject', 0);
    end
    uncuedReward = uncuedTrials & rewardTrials;
    uncuedPunish = uncuedTrials & punishTrials;
    uncuedShock = uncuedTrials & shockTrials;
    trialTypes = 1:max(TE.trialType);
    trialsByType = cell(size(trialTypes));
    for ttc = 1:max(TE.trialType)
        trialsByType{ttc} = filterTE(TE, 'trialType', trialTypes(ttc), 'reject', 0);
    end
    trialCount = [1:length(TE.filename)]';
    clear ttc;
