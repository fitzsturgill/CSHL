% ChAT_26 + pupilometry analysis script

% exclude: Aug20- bad pupil, Aug 21, no photomtery

% sessions = bpLoadSessions;
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug16_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug17_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug18_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug19_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug22_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug23_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***
% ***Loaded ChAT_26_CuedOutcome_odor_complete_Aug24_2016_Session1.mat from E:\Data\ChAT_26\CuedOutcome_odor_complete\Session Data\***

TE = makeTE_CuedOutcome_Odor_Complete(sessions);
TE.Photometry = processTrialAnalysis_Photometry2(sessions);
TE = addPupilometryToTE(TE);
%%
% TE.pupBaseline = bpCalcPeak_Pupil(TE.pupil, [-5 -3], TE.Us);
% TE.pupPeak = bpCalcPeak_Pupil(TE.pupil, [1 1.5], TE.Us);
% TE.phPeak = bpCalcPeak_dFF(TE.Photometry, 1, [0.1 0.4], TE.Us);
% TE.csPeak = bpCalcPeak_dFF(TE.Photometry, 1, [-1 0], TE.Us);

TE.pupBaseline = bpCalcPeak_Pupil(TE.pupil, [-5 -3] + 6, []);
TE.pupPeak = bpCalcPeak_Pupil(TE.pupil, [1 1.5] + 6, []);
TE.phPeak = bpCalcPeak_dFF(TE.Photometry, 1, [0.1 0.4] + 6, []);
TE.csPeak = bpCalcPeak_dFF(TE.Photometry, 1, [-2 -1] + 6, []);
TE.pupCsPeak = bpCalcPeak_Pupil(TE.pupil, [-1 0] + 6, []);


%% use TrialCutoffs to add reject field to TE
% TrialCutoffs = [110 85 90 140 100 80 95 135]; %
TrialCutoffs = [110 85 90 140 80 95 135]; %
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
    includedSessions = 1:7;
    highValueTrials = filterTE(TE, 'trialType', 1:3, 'reject', 0, 'sessionIndex', includedSessions);
    lowValueTrials = filterTE(TE, 'trialType', 4:6, 'reject', 0, 'sessionIndex', includedSessions);
    rewardTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0, 'sessionIndex', includedSessions);
    punishTrials = filterTE(TE, 'trialOutcome', 2, 'reject', 0, 'sessionIndex', includedSessions);
    omitTrials = filterTE(TE, 'trialOutcome', 3, 'reject', 0, 'sessionIndex', includedSessions);
    uncued = filterTE(TE, 'trialType', 7:9, 'reject', 0, 'sessionIndex', includedSessions);
    trialTypes = 1:9;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0, 'sessionIndex', includedSessions);
    end
    allTrials = filterTE(TE, 'reject', 0);
%%
ensureFigure('test', 1);
subplot(1,2,1);
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{1}, :)), 'c'); % hi reward
hold on;
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{7}, :)), 'b'); % uncued reward
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{5}, :)), 'r'); % lo punish
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{9}, :)), 'k'); % uncued neutral
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(allTrials, :)), 'y'); % uncued neutral
set(gca, 'YLim', [0.7 1.3], 'XGrid', 'on');

subplot(1,2,2);
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{1}, :)), 'b'); % hi reward
hold on;
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{2}, :)), 'r'); % hi punish
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{3}, :)), 'k'); % hi omission
set(gca, 'YLim', [0.7 1.3]);

%%
ensureFigure('averages_by_outcome', 1);
subplot(2,2,1);
phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1, 'linespec', {'c', 'b', 'k'});
set(gca, 'XLim', [-4 4], 'XGrid', 'on');set(gca, 'FontSize', 12, 'TickDir', 'out');
title('ChAT dF/F');
ylabel('dFF');
subplot(2,2,2); 
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{1}, :)), 'c'); % 
hold on;
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{4}, :)), 'b'); % 
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{7}, :)), 'k'); % 
set(gca, 'YLim', [0.7 1.3], 'XLim', [-4 4], 'XGrid', 'on');set(gca, 'FontSize', 12, 'TickDir', 'out');
title('Pupil');
ylabel('Pup Diam. norm');

subplot(2,2,3);
phPlotAverageFromTE(TE, trialsByType([5 2 8]), 1, 'linespec', {'m', 'r', 'k'});
set(gca, 'XLim', [-4 4], 'XGrid', 'on');set(gca, 'FontSize', 12, 'TickDir', 'out');
xlabel('time from reinforcement (s)');
ylabel('dFF');
subplot(2,2,4); 
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{5}, :)), 'm'); % 
hold on;
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{2}, :)), 'r'); % 
plot(TE.pupil.xData, nanmean(TE.pupil.pupDiameterNorm(trialsByType{8}, :)), 'k'); % 
set(gca, 'YLim', [0.7 1.3], 'XLim', [-4 4], 'XGrid', 'on');set(gca, 'FontSize', 12, 'TickDir', 'out');
xlabel('time from reinforcement (s)');
ylabel('Pup Diam. norm');

%%

pupilZero = bpX2pnt(6, 60, 0);
% phZero = bpX2pnt(6, 20, 0);
ensureFigure('averages_by_cue', 1);
subplot(1,2,1);
phPlotAverageFromTE(TE, {highValueTrials, lowValueTrials}, 1, 'linespec', {'b', 'r'}, 'window', [-4 0]); hold on;
set(gca, 'XLim', [-4 0], 'XGrid', 'on');set(gca, 'FontSize', 12, 'TickDir', 'out');
title('ChAT dF/F (hi vs lo value cue)');
subplot(1,2,2); 
plot(TE.pupil.xData(1:pupilZero), nanmean(TE.pupil.pupDiameterNorm(highValueTrials, (1:pupilZero))), 'b'); % 
hold on;
plot(TE.pupil.xData(1:pupilZero), nanmean(TE.pupil.pupDiameterNorm(lowValueTrials, (1:pupilZero))), 'r'); % 
set(gca, 'YLim', [0.9 1.15], 'XLim', [-4 0], 'XGrid', 'on');set(gca, 'FontSize', 12, 'TickDir', 'out');
xlabel('time from reinforcement (s)');
ylabel('Pupil Diameter norm (by session)');
title('Pupil (hi vs lo value cue)');

%% 
ensureFigure('pupBaseline_vs_phPeak', 1);
for counter = 1:9
    subplot(3,3,counter);
    scatter(TE.pupBaseline.data(trialsByType{counter}), TE.phPeak.data(trialsByType{counter}));
end

%% 
ensureFigure('cue_vs_cue', 1);
scatter(TE.pupCsPeak.data(highValueTrials), TE.csPeak.data(highValueTrials), 'MarkerFaceColor', 'r'); hold on;
scatter(TE.pupCsPeak.data(uncued), TE.csPeak.data(uncued), 'MarkerFaceColor', 'k'); 
xlabel('Pupil Diameter (cue)'); ylabel('dFF (cue)');
title('high value cue vs dummy cue');
set(gca, 'FontSize', 12, 'TickDir', 'out');







