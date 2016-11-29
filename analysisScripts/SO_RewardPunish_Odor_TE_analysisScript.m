% SO_RewardPunish_Odor_TE_analysisScript

saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_SO_RewardPunish_Odor_v2(sessions);
%%
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial', 'zeroField', 'PostUsRecording');

%%
basepath = uigetdir;
subjectName = TE.filename{1}(1:7);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

%%
TE.csLicks = countEventFromTE(TE, 'Port1In', [-2 0], TE.PostUsRecording);
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.PostUsRecording);
truncateSessionsFromTE(TE, 'init');
%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% generate trial lookups for different combinations of conditions
    validTrials = filterTE(TE, 'reject', 0);
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

    
    
%%
h=ensureFigure('Averages', 1); 
pm = [1 1];
subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 5]), 1,...
    'linespec', {'c', 'b'}, 'window', [-0.25 2]); % cued, unexpected reward
legend(hl, {'cued', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
% title('Photometry'); ylabel('dF/F'); xlabel('time from reinforcement (s)');
textBox(TE.filename{1}(1:7));
set(gca, 'XLim', [-0.25 2], 'TickDir', 'out');
%%
if saveOn
    saveas(gcf, fullfile(savepath, 'ChAT_13_averages.fig'));
    saveas(gcf, fullfile(savepath, 'ChAT_13_averages.jpg'));
    saveas(gcf, fullfile(savepath, 'ChAT_13_averages.epsc'));    
end
