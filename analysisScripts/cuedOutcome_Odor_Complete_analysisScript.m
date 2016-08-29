
% using bpLoadSessions loaded the following data:
% ChAT_32, August 15 -> 26
% ChAT_31, August 15 -> 24
% ChAT_30, August 17 -> 24
% ChAT_29, August 17 -> 24
% ChAT_28, August 17 -> 23
% ChAT_26, August 17 -> 24
% and used makeTE_CuedOutcome_Odor_Complete to make a TE structure where each 
% element in the structure was a TE corresponding to a single mouse

savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
%% Make tiled array of licks in receipt of reward to visualize satiation/ lapsing behavior towards end of each session
ensureFigure('RewardLickRate_crossSessions', 1);
for tei = 1:6
    rewardTrials = filterTE(TE(tei), 'trialOutcome', 1);
    subplot(3,2,tei); plot(find(rewardTrials), smooth(TE(tei).usLicks.rate(rewardTrials), 5)); hold on; 
    plot(find(rewardTrials), [0; diff(TE(tei).sessionIndex(rewardTrials))] * 10);
    ylabel('Lick/s, Us'); xlabel('trial #'); textBox(TE(tei).filename{1}(1:7));
    set(gca, 'YLim', [0 10]);
end
saveas(gcf, fullfile(savepath, 'RewardLickRate_crossSessions.fig'));

%% make tiled array of antic. licks for low and high value odors vs trial number
smoothFactor = 11;
ensureFigure('RewardLickRate_crossSessions', 1);
for tei = 1:6
    highTrials = filterTE(TE(tei), 'trialType', 1:3);
    lowTrials = filterTE(TE(tei), 'trialType', 4:6);
    rewardTrials = filterTE(TE(tei), 'trialOutcome', 1);    
    subplot(3,2,tei); plot(find(highTrials), smooth(TE(tei).csLicks.rate(highTrials), smoothFactor), 'b'); hold on; 
    plot(find(lowTrials), smooth(TE(tei).csLicks.rate(lowTrials), smoothFactor), 'r');
    plot(find(rewardTrials), smooth(TE(tei).usLicks.rate(rewardTrials), smoothFactor), 'k')    
    plot(1:length(TE(tei).trialNumber), [0; diff(TE(tei).sessionIndex)] * 10, 'g');
    ylabel('Lick/s, Us'); xlabel('trial #'); textBox(TE(tei).filename{1}(1:7));
    set(gca, 'YLim', [0 10]);
end
saveas(gcf, fullfile(savepath, 'AnticLickRate_crossSessions.fig'));



%%

%     highValueTrials = filterTE(TE(tei), 'trialType', 1:3);
%     lowValueTrials = filterTE(TE(tei), 'trialType', 4:6);
%     rewardTrials = filterTE(TE(tei), 'trialOutcome', 1);
%     punishTrials = filterTE(TE(tei), 'trialOutcome', 2);
%     omitTrials = filterTE(TE(tei), 'trialOutcome', 3);
%     ensureFigure('test', 1); plot(find(rewardTrials), smooth(TE(tei).usLicks.rate(rewardTrials), 5)); hold on; 
%     plot(find(rewardTrials), [0; diff(TE(tei).sessionIndex(rewardTrials))] * 10);