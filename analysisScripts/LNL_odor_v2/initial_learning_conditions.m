LNL_conditions;
% find first reversal if it occurs on day 2+
day1Trials = TE.sessionIndex == 1 & validTrials;
day2PlusTrials = TE.sessionIndex >= 2 & validTrials;
firstReversalBlock = TE.BlockNumber(find(TE.sessionIndex >= 2 & TE.OdorValveIndex == 1 & TE.CSValence == -1, 1)); % block number where 1st odor has a negative valence.
if ~isempty(firstReversalBlock)
    firstReversalTrial = find(TE.BlockNumber == firstReversalBlock, 1);
    day2PlusTrials(firstReversalTrial:end) = false; % if firstReversalTrial is empty nothing happens...
end