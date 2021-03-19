% script to generate trial lookups by condition for
% cuedOutcome_Odor_complete

%% generate trial lookups for different combinations of conditions
%     validTrials = filterTE(TE, 'reject', 0);
%     highValueTrials = filterTE(TE, 'trialType', 1:3, 'reject', 0);
%     lowValueTrials = filterTE(TE, 'trialType', 4:6, 'reject', 0);
%     uncuedTrials = filterTE(TE, 'trialType', 7:9, 'reject', 0);    
%     rewardTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
%     punishTrials = filterTE(TE, 'trialOutcome', 2, 'reject', 0);    
%     omitTrials = filterTE(TE, 'trialOutcome', 3, 'reject', 0);
    validTrials = filterTE(TE, 'reject', 0);
%     badTrials1 = cellfun(@(x,y) y(end) - x(1), TE.PreCsRecording, TE.foreperiod) < 3.8;
    badTrials1 = cellfun(@(x,y) y(end) - x(1), TE.PreCsRecording, TE.foreperiod) ~= 4;
    badTrials2 = cellfun(@(x) x(end) - x(1), TE.Cue) > 1;
    badTrials3 = cellfun(@(x) x(end) - x(1), TE.PostUsRecording) > 6;
    badTrials4 = cellfun(@(x) x(end) - x(1), TE.Reward) > 0.3; 
    badTrials5 = cellfun(@(x) x(end) - x(1), TE.Punish) > 0.6;
    badTrials6 = cellfun(@(x) x(end) - x(1), TE.Neutral) > 0.3;
%     badTrials7 = TE.RT < 0.05;
    badTrials = badTrials1 + badTrials2 + badTrials3 + badTrials4  + badTrials5  + badTrials6;
    allTrials = filterTE(TE, 'reject', 0) & ~badTrials;
    fpLickTrials = TE.fpLicks.count > 0;
%     badTrials = badTrials1 + badTrials2 + badTrials3 + badTrials4  + badTrials5  + badTrials6 + badTrials7 + fpLickTrials;
    
    Sound1Trials = filterTE(TE, 'SoundValveIndex', 1, 'reject', 0) & ~badTrials;
    Sound2Trials = filterTE(TE, 'SoundValveIndex', 2, 'reject', 0) & ~badTrials; 
    Sound3Trials = filterTE(TE, 'SoundValveIndex', 3, 'reject', 0) & ~badTrials;
%     Sound4Trials = filterTE(TE, 'SoundValveIndex', 4, 'reject', 0) & ~badTrials;
    uncuedTrials = filterTE(TE, 'SoundValveIndex', 0, 'reject', 0) & ~badTrials;
%     anticipTrials = TE.csLicks.count >= 2;
%     noanticipTrials = TE.csLicks.count < 2;     

%     rewardTrials_SL = filterTE(TE, 'ReinforcementOutcome', 'Reward', 'reject', 0) & ~badTrials;
%     punishTrials_SL = filterTE(TE, 'ReinforcementOutcome', 'Punish', 'reject', 0) & ~badTrials;
%     omissionTrials_SL = filterTE(TE, 'ReinforcementOutcome', 'Neutral', 'reject', 0) & ~badTrials;         
 
%     Sound1Reward = filterTE(TE, 'SoundValveIndex', 1, 'ReinforcementOutcome', 'Reward', 'reject', 0) & ~badTrials;
%     Sound1Omission = filterTE(TE, 'SoundValveIndex', 1, 'ReinforcementOutcome', 'Neutral', 'reject', 0) & ~badTrials;
%     Sound2Reward = filterTE(TE, 'SoundValveIndex', 2, 'ReinforcementOutcome', 'Reward', 'reject', 0) & ~badTrials;
%     Sound2Punish = filterTE(TE, 'SoundValveIndex', 2, 'ReinforcementOutcome', 'Punish', 'reject', 0) & ~badTrials;
%     Sound2Omission = filterTE(TE, 'SoundValveIndex', 2, 'ReinforcementOutcome', 'Neutral', 'reject', 0) & ~badTrials;
%     Sound3Reward = filterTE(TE, 'SoundValveIndex', 3, 'ReinforcementOutcome', 'Reward', 'reject', 0) & ~badTrials;
%     Sound3Omission = filterTE(TE, 'SoundValveIndex', 3, 'ReinforcementOutcome', 'Neutral', 'reject', 0) & ~badTrials;
    uncuedReward = filterTE(TE, 'SoundValveIndex', 0, 'ReinforcementOutcome', 'Reward', 'reject', 0) & ~badTrials;
    uncuedPunish = filterTE(TE, 'SoundValveIndex', 0, 'ReinforcementOutcome', 'Punish', 'reject', 0) & ~badTrials;

    hitTrials = filterTE(TE, 'SoundValveIndex', 1, 'LickAction', 'lick', 'reject', 0) & ~badTrials;
    missTrials = filterTE(TE, 'SoundValveIndex', 1, 'LickAction', 'nolick', 'reject', 0) & ~badTrials;
    FATrials = filterTE(TE, 'SoundValveIndex', 2, 'LickAction', 'lick', 'reject', 0) & ~badTrials;
    CRTrials = filterTE(TE, 'SoundValveIndex', 2, 'LickAction', 'nolick', 'reject', 0) & ~badTrials;
    
    Sound1_50_Trials = filterTE(TE, 'SoundValveIndex', 1, 'SoundAmplitude', 50, 'reject', 0) & ~badTrials;
    Sound1_40_Trials = filterTE(TE, 'SoundValveIndex', 1, 'SoundAmplitude', 40, 'reject', 0) & ~badTrials;
    Sound1_30_Trials = filterTE(TE, 'SoundValveIndex', 1, 'SoundAmplitude', 30, 'reject', 0) & ~badTrials;
    Sound1_20_Trials = filterTE(TE, 'SoundValveIndex', 1, 'SoundAmplitude', 20, 'reject', 0) & ~badTrials;
    Sound1_10_Trials = filterTE(TE, 'SoundValveIndex', 1, 'SoundAmplitude', 10, 'reject', 0) & ~badTrials;
    Sound2_50_Trials = filterTE(TE, 'SoundValveIndex', 2, 'SoundAmplitude', 50, 'reject', 0) & ~badTrials;
    Sound2_40_Trials = filterTE(TE, 'SoundValveIndex', 2, 'SoundAmplitude', 40, 'reject', 0) & ~badTrials;
    Sound2_30_Trials = filterTE(TE, 'SoundValveIndex', 2, 'SoundAmplitude', 30, 'reject', 0) & ~badTrials;
    Sound2_20_Trials = filterTE(TE, 'SoundValveIndex', 2, 'SoundAmplitude', 20, 'reject', 0) & ~badTrials;
    Sound2_10_Trials = filterTE(TE, 'SoundValveIndex', 2, 'SoundAmplitude', 10, 'reject', 0) & ~badTrials;

    hit50Trials = filterTE(TE, 'SoundValveIndex', 1, 'LickAction', 'lick', 'SoundAmplitude', 50, 'reject', 0) & ~badTrials;
    hit40Trials = filterTE(TE, 'SoundValveIndex', 1, 'LickAction', 'lick', 'SoundAmplitude', 40, 'reject', 0) & ~badTrials;
    hit30Trials = filterTE(TE, 'SoundValveIndex', 1, 'LickAction', 'lick', 'SoundAmplitude', 30, 'reject', 0) & ~badTrials;
    hit20Trials = filterTE(TE, 'SoundValveIndex', 1, 'LickAction', 'lick', 'SoundAmplitude', 20, 'reject', 0) & ~badTrials;
    hit10Trials = filterTE(TE, 'SoundValveIndex', 1, 'LickAction', 'lick', 'SoundAmplitude', 10, 'reject', 0) & ~badTrials;
% 
    FA50Trials = filterTE(TE, 'SoundValveIndex', 2, 'LickAction', 'lick', 'SoundAmplitude', 50, 'reject', 0) & ~badTrials;
    FA40Trials = filterTE(TE, 'SoundValveIndex', 2, 'LickAction', 'lick', 'SoundAmplitude', 40, 'reject', 0) & ~badTrials;
    FA30Trials = filterTE(TE, 'SoundValveIndex', 2, 'LickAction', 'lick', 'SoundAmplitude', 30, 'reject', 0) & ~badTrials;
    FA20Trials = filterTE(TE, 'SoundValveIndex', 2, 'LickAction', 'lick', 'SoundAmplitude', 20, 'reject', 0) & ~badTrials;
    FA10Trials = filterTE(TE, 'SoundValveIndex', 2, 'LickAction', 'lick', 'SoundAmplitude', 10, 'reject', 0) & ~badTrials;
%     hit50Trials = filterTE(TE, 'trialType', [1], 'LickAction', 'lick', 'reject', 0) & ~badTrials;
%     hit40Trials = filterTE(TE, 'trialType', [2], 'LickAction', 'lick', 'reject', 0) & ~badTrials;
%     hit30Trials = filterTE(TE, 'trialType', [3], 'LickAction', 'lick', 'reject', 0) & ~badTrials;
% 
%     FA50Trials = filterTE(TE, 'trialType', [4], 'LickAction', 'lick', 'reject', 0) & ~badTrials;
%     FA40Trials = filterTE(TE, 'trialType', [5], 'LickAction', 'lick', 'reject', 0) & ~badTrials;
%     FA30Trials = filterTE(TE, 'trialType', [6], 'LickAction', 'lick', 'reject', 0) & ~badTrials;
%     
    Type1Trials = filterTE(TE, 'trialType', [1], 'reject', 0) & ~badTrials;
    Type2Trials = filterTE(TE, 'trialType', [2], 'reject', 0) & ~badTrials;
    Type3Trials = filterTE(TE, 'trialType', [3], 'reject', 0) & ~badTrials;
    Type4Trials = filterTE(TE, 'trialType', [4], 'reject', 0) & ~badTrials;
    Type5Trials = filterTE(TE, 'trialType', [5], 'reject', 0) & ~badTrials;
    Type6Trials = filterTE(TE, 'trialType', [6], 'reject', 0) & ~badTrials;
    Type7Trials = filterTE(TE, 'trialType', [7], 'reject', 0) & ~badTrials;
    Type8Trials = filterTE(TE, 'trialType', [8], 'reject', 0) & ~badTrials;
    Type9Trials = filterTE(TE, 'trialType', [9], 'reject', 0) & ~badTrials; 
%     
%     uncuedReward = Type6Trials;
%     uncuedPunish = Type7Trials;
    
    firstNTrials = 10;
    
    nSessions = max(TE.sessionIndex);
    Sound1FirstNTrials = [];
    Sound2FirstNTrials = [];
    for ttc = 1:nSessions
        trialsThisSession = find(Sound1Trials & (TE.sessionIndex == ttc));
        Sound1FirstNTrials = [Sound1FirstNTrials; trialsThisSession(1:min(firstNTrials, length(trialsThisSession)))];
        trialsThisSession = find(Sound2Trials & (TE.sessionIndex == ttc));
        Sound2FirstNTrials = [Sound2FirstNTrials; trialsThisSession(1:min(firstNTrials, length(trialsThisSession)))];        
    end
        
    trialTypes = 1:max(TE.trialType);
    trialsByType = cell(size(trialTypes));
    for ttc = 1:length(trialTypes)
        trialsByType{ttc} = filterTE(TE, 'trialType', trialTypes(ttc), 'reject', 0);
    end

    trialCount = [1:length(TE.filename)];
    
    clear ttc;
