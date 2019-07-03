% FrankenLNL_varyRewardSize_conditions

% script to generate trial lookups by condition for
% frankenLNL_4odors, block 10, rewardPunishBlocks

%% generate trial lookups for different combinations of conditions

validTrials = filterTE(TE, 'reject', 0);
largeRewardTrials = filterTE(TE, 'trialType', 1, 'reject', 0);
mediumRewardTrials = filterTE(TE, 'trialType', 2, 'reject', 0);
smallRewardTrials = filterTE(TE, 'trialType', 3, 'reject', 0);
neutralTrials = filterTE(TE, 'trialType', 4, 'reject', 0);

trialTypes = 1:max(TE.trialType);
trialsByType = cell(size(trialTypes));
for ttc = 1:max(TE.trialType)
    trialsByType{ttc} = filterTE(TE, 'trialType', trialTypes(ttc), 'reject', 0);
end
trialCount = [1:length(TE.filename)]';
clear ttc;
