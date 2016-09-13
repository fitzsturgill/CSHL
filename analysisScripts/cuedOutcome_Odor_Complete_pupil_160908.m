% ChAT_26 + pupilometry analysis script


% sessions = bpLoadSessions;
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug16_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug17_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug18_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug19_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug20_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug22_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug23_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug24_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***

TE = makeTE_CuedOutcome_Odor_Complete(sessions);
TE.Photometry = processTrialAnalysis_Photometry2(sessions);
TE = addPupilometryToTE(TE);


%% use TrialCutoffs to add reject field to TE
TrialCutoffs = [110 85 90 140 100 80 95 135]; %
TE = ensureField(TE, 'reject', 'mat');

nSessions = length(unique(TE.sessionIndex));
if nSessions == length(TrialCutoffs)
    reject = ones(size(TE.trialNumber));
    for i = 1:nSessions
        include = (TE.sessionIndex == i) & (TE.trialNumber' <= TrialCutoffs(i));
        reject(include) = 0;
    end
    TE.reject = reject;
end

%% generate trial lookups for different combinations of conditions
    includedSessions = 5;
    highValueTrials = filterTE(TE, 'trialType', 1:3, 'reject', 0, 'sessionIndex', includedSessions);
    lowValueTrials = filterTE(TE, 'trialType', 4:6, 'reject', 0, 'sessionIndex', includedSessions);
    rewardTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0, 'sessionIndex', includedSessions);
    punishTrials = filterTE(TE, 'trialOutcome', 2, 'reject', 0, 'sessionIndex', includedSessions);
    omitTrials = filterTE(TE, 'trialOutcome', 3, 'reject', 0, 'sessionIndex', includedSessions);
    trialTypes = 1:9;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0, 'sessionIndex', includedSessions);
    end
%%
ensureFigure('test', 1);
subplot(1,2,1);
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{1}, :)), 'c'); % hi reward
hold on;
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{7}, :)), 'b'); % uncued reward
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{5}, :)), 'r'); % lo punish
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{9}, :)), 'k'); % uncued neutral
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(:, :)), 'y'); % uncued neutral
set(gca, 'YLim', [0.5 1.5]);

subplot(1,2,2);
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{1}, :)), 'b'); % hi reward
hold on;
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{2}, :)), 'r'); % hi punish
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{3}, :)), 'k'); % hi omission
set(gca, 'YLim', [0.5 1.5]);





