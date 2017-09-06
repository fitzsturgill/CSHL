% CCN Talk Script

%% desktop
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\';
saveOn = 1;

%% laptop
savepath = 'C:\Users\Fitz\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk';
saveOn = 1;
%% Snippet to make lick histogram for graded value task from DC_35: (or DC_26?)

 % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-6 4], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};

    ensureFigure('Lick_Hist_CCN', 1); axes;

    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 4 7]), 'Port1In', varargin{:});
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    set(gca, 'XLim', [-6 4], 'YLim', [-2 15]);
    addStimulusPatch(gca, [-3 -2 10], '', [0.7 0.7 0.7]) 
    addStimulusPatch(gca, [-0.1 0.1], '');
    ylabel('Licks (Hz)'); xlabel('Time from reinforcement (s)');
    formatFigureTalk([4 3]);
    if saveOn    
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN.fig'));
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN.jpg'));    
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN.meta'));                   
    end
    
%% Us and Cs scatter plots
CuedOutcome_pooledAnalysis_script2;
  
%% Snippet to make dopamine photometry averages from DC_26:

    ensureFigure('Ph_Hist_Dopamine_CCN', 1); axes('FontSize', 12, 'LineWidth', 1); 
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1, 'FluorDataField', 'ZS', 'linespec', {'b', 'r', 'g'}); % reward, varying degrees of expectation
    set(gca, 'XLim', [-4 4]); ylabel('Z Score'); xlabel('time from reinforcement (s)');     

    if saveOn    
        saveas(gcf, fullfile(savepath, 'Ph_Hist_Dopamine_CCN.fig'));
        saveas(gcf, fullfile(savepath, 'Ph_Hist_Dopamine_CCN.jpg'));    
        saveas(gcf, fullfile(savepath, 'Ph_Hist_Dopamine_CCN.meta'));           
    end
    
    
%% Snippet to make cuedOutcome cue response photometry histogram from ChAT_42

    ensureFigure('CuedOutcome_Cue', 1); 
    [ha, hl] = phPlotAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 1,...
        'window', [-4 3], 'linespec', {'b', 'r', 'g'}, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue);
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    set(gca, 'XLim', [-4 3]);
    addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7]) 
    ylabel('Fluorescence (Z Score)'); xlabel('Time from cue (s)');
    formatFigureTalk([4 3]);    

    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue.fig'));
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue.jpg'));    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue.meta'));           
    end
    
%% Snippet to make cuedOutcome cue response Lick histogram from ChAT_42
 % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-4 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};
    ensureFigure('CuedOutcome_Cue_Lick', 1);
    [ha, hl] = plotEventAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 'Port1In', varargin{:});
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    set(gca, 'XLim', [-4 3]);
    addStimulusPatch(gca, [0 1], 'odor', [0.7 0.7 0.7]) 
    ylabel('Licks (Hz)'); xlabel('Time from cue (s)');
    formatFigureTalk([4 2]);        
    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_Lick.fig'));
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_Lick.jpg'));    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_Lick.meta'));           
    end    
    
%% Snippet to make complete photometry histogram (reward condition from ChAT_42    
    
    ensureFigure('CuedOutcome_Reward', 1); 
    [ha, hl] = phPlotAverageFromTE(TE, {trialsByType{1}, trialsByType{4}, trialsByType{7}}, 1,...
        'window', [-2 2], 'linespec', {'b', 'r', 'g'}, 'FluorDataField', 'ZS');
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    set(gca, 'XLim', [-2 2], 'YLim', [-1.2 2.5]);
%     addStimulusPatch(gca, [-3 -2 10], '', [0.7 0.7 0.7]) 
    addStimulusPatch(gca, [-0.1 0.1], '');
    ylabel('Fluorescence (Z Score)'); xlabel('Time from reinforcement (s)');
    formatFigureTalk([4 3]);    

    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward.fig'));
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward.jpg'));    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward.meta'));           
    end
    
    %% Snippet to make cuedOutcome Reward response Lick histogram from ChAT_42
 % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-2 2], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};
    ensureFigure('CuedOutcome_Reward_Lick', 1);
    [ha, hl] = plotEventAverageFromTE(TE, {trialsByType{1}, trialsByType{4}, trialsByType{7}}, 'Port1In', varargin{:});
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    set(gca, 'XLim', [-2 2], 'Box', 'off');
    addStimulusPatch(gca, [-0.1 0.1], '');
    ylabel('Licks (Hz)'); xlabel('Time from reinforcement (s)');
    formatFigureTalk([4 2]);        
    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward_Lick.fig'));
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward_Lick.jpg'));    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward_Lick.meta'));           
    end    
    
    
    
%%
%% Lick and Ph rasters from DC_20- by CS valence
CLimFactor = 2;
CLimFactor2 = 3;
trialRange = [170 340]; % global trial range
% trialRange = [560 680]; % global trial range
includedTrials = zeros(size(csPlusTrials));
includedTrials(trialRange(1):trialRange(2)) = 1;
includedTrials = logical(includedTrials);
reversals = find(diff(TE.BlockNumber(csPlusTrials & includedTrials & rewardTrials)));

    saveName = ['reversals_phRasters_byValence' num2str(trialRange(1)) '_to_' num2str(trialRange(2))];
    h=ensureFigure(saveName, 1);
%     mcPortraitFigSetup(h);
    

%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,3,1); 
    eventRasterFromTE(TE, csPlusTrials & includedTrials & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
%     title('CS+'); ylabel('trial number');
    set(gca, 'XLim', [-4 7], 'TickDir', 'Out'); 
    set(gca, 'FontSize', 10)
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'k', 'LineWidth', 2); % reversal lines    
    set(gca, 'XLim', [-2 6]); ylabel('Trial #');
    addStimulusBar(gca, [0 1 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [2.9 3.1 0], '', [0.3 0.3 0.3], 5);
    
    subplot(1,3,2); phRasterFromTE(TE, csPlusTrials & includedTrials & rewardTrials, 1, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 10, 'TickDir', 'Out')
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % reversal lines        
        set(gca, 'XLim', [-2 6], 'YTick', []); xlabel('Time from cue (s)');
    addStimulusBar(gca, [0 1 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [2.9 3.1 0], '', [0.3 0.3 0.3], 5); 
     
    subplot(1,3,3); phRasterFromTE(TE, csPlusTrials & includedTrials & rewardTrials, 2, 'CLimFactor', CLimFactor2, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 10, 'TickDir', 'Out')
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % reversal lines       
        set(gca, 'XLim', [-2 6]);
    set(gca, 'XLim', [-2 6], 'YTick', []);
    addStimulusBar(gca, [0 1 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [2.9 3.1 0], '', [0.3 0.3 0.3], 5);    
    formatFigureTalk([4.5 3]);        
  


    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.meta']));   
    end

    %% Lick and Ph rasters from DC_20- by Odor
CLimFactor = 2;
CLimFactor2 = 3;
trialRange = [170 340]; odorTrials = Odor2Trials; % global trial range
% trialRange = [568 680]; odorTrials = Odor2Trials; % global trial range

includedTrials = zeros(size(csPlusTrials));
includedTrials(trialRange(1):trialRange(2)) = 1;
includedTrials = logical(includedTrials);
reversals = find(diff(TE.BlockNumber(odorTrials & includedTrials)));

    saveName = ['reversals_phRasters_byOdor' num2str(trialRange(1)) '_to_' num2str(trialRange(2))];
    h=ensureFigure(saveName, 1);
%     mcPortraitFigSetup(h);
    

%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,3,1); 
    eventRasterFromTE(TE, odorTrials & includedTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
%     title('CS+'); ylabel('trial number');
    set(gca, 'FontSize', 10)
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    set(gca, 'XLim', [-2 6]); ylabel('Trial #');
    addStimulusBar(gca, [0 1 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [2.9 3.1 0], '', [0.3 0.3 0.3], 5);
    
    subplot(1,3,2); phRasterFromTE(TE, odorTrials & includedTrials, 1, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines        
        set(gca, 'XLim', [-2 6], 'YTick', []); xlabel('Time from cue (s)');
    addStimulusBar(gca, [0 1 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [2.9 3.1 0], '', [0.3 0.3 0.3], 5);        
     
    subplot(1,3,3); phRasterFromTE(TE, odorTrials & includedTrials, 2, 'CLimFactor', CLimFactor2, 'trialNumbering', 'consecutive');
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines       
    set(gca, 'XLim', [-2 6], 'YTick', []);
    addStimulusBar(gca, [0 1 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [2.9 3.1 0], '', [0.3 0.3 0.3], 5);    
    formatFigureTalk([4.5 3]);


    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.meta']));   
    end
    
%% sorting by behavior (doesn't really work well)


excludeFirstNTrials = 20;
excludeFirstTrials = TE.trialNumber > excludeFirstNTrials;
threshPercentile = 0.25;

TE.csLicksAll = countEventFromTE(TE, 'Port1In', [0 3], TE.Cue);
ensureFigure('cueLick_Hist', 1);
[sorted index] = cum(TE.csLicksAll.rate(csPlusTrials & excludeFirstTrials));
plot(sorted, index); xlabel('lick rate');
threshRate = percentile(TE.csLicksAll.rate(csPlusTrials & excludeFirstTrials), threshPercentile);
hitTrialsRewarded = csPlusTrials & excludeFirstTrials & TE.csLicksAll.rate > threshRate & rewardTrials;
missTrialsRewarded = csPlusTrials & excludeFirstTrials & TE.csLicksAll.rate <= threshRate & rewardTrials;

ensureFigure('rev_phAveragesByBehavior', 1);
subplot(1, 2, 1);
[ha, hl] = phPlotAverageFromTE(TE, {hitTrialsRewarded, missTrialsRewarded}, 1,...
'FluorDataField', 'ZS', 'window', [-1, 5], 'linespec', {'g', 'r'}); %high value, reward
legend(hl, {'hit', 'miss'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
ylabel('CS+, outcomes'); set(gca, 'XLim', [-1, 5]); title('ACh.');%set(gca, 'YLim', ylim);

subplot(1,2,2);
[ha, hl] = phPlotAverageFromTE(TE, {hitTrialsRewarded, missTrialsRewarded}, 2,...
'FluorDataField', 'ZS', 'window', [-1, 5], 'linespec', {'g', 'r'}); %high value, reward
legend(hl, {'hit', 'miss'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
set(gca, 'XLim', [-1, 5]); title('Dop.')%set(gca, 'YLim', ylim);
formatFigureTalk([6,3]);

%% sorting by number of cue reward pairings post reversal

excludeFirstNTrials = 20;
excludeFirstTrials = TE.trialNumber > excludeFirstNTrials;
rewardPairingsThresh = 10; % 
nSessions = max(TE.sessionIndex);
trialMatrix = zeros(length(TE.filename), nSessions);
for counter = 1:nSessions
    trialsThisSession = TE.sessionIndex == counter;
    thisReverse = find(TE.BlockChange & trialsThisSession);
    postReversalTrials = trialsThisSession; postReversalTrials(1:thisReverse - 1) = 0;
    pairingsThisReverse = postReversalTrials & csPlusTrials & rewardTrials;    
    trialSearch = find(pairingsThisReverse);
    threshTrial = trialSearch(rewardPairingsThresh);
    firstTrials = postReversalTrials; firstTrials(threshTrial:end) = 0;
    trialMatrix(:,counter) = firstTrials;
end
periReversalTrials = any(trialMatrix, 2);
extraReversalTrials = ~periReversalTrials;

patchHue = [0.6 0.6 0.6];
saveName = 'rev_phAveragesByPairings_ChAT';
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & periReversalTrials, csPlusTrials & rewardTrials & extraReversalTrials}, 1,...
'FluorDataField', 'ZS', 'window', [-1, 5], 'linespec', {'g', 'b'}); %high value, reward
legend(hl, {['\bf\color{green} Initial trials'], '\bf\color{blue} Others'}, 'Location', 'northwest', 'FontSize', 16, 'Box', 'on'); 
ylabel('Fluor. (Zscore)'); set(gca, 'XLim', [-1, 5], 'YLim', [-.5 3]);
xlabel('Time from odor (s)');
addStimulusPatch(gca, [0 1], '', patchHue);
addStimulusPatch(gca, [2.9 3.1], '', patchHue);
formatFigureTalk;
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.meta']));   
end
    

saveName = 'rev_phAveragesByPairings_DAT';
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & periReversalTrials, csPlusTrials & rewardTrials & extraReversalTrials}, 2,...
'FluorDataField', 'ZS', 'window', [-1, 5], 'linespec', {'g', 'b'}); %high value, reward
% legend(hl, {['\bf\color{green} Initial'], '\bf\color{blue} Others'}, 'Location', 'northwest', 'FontSize', 16, 'Box', 'off'); 
set(gca, 'XLim', [-1, 5], 'YLim', [-.5 3]);
ylabel('Fluor. (Zscore)');
xlabel('Time from odor (s)');
addStimulusPatch(gca, [0 1], '', patchHue);
addStimulusPatch(gca, [2.9 3.1], '', patchHue);
formatFigureTalk;

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.meta']));   
end

%% as above, outcome only

patchHue = [0.6 0.6 0.6];
saveName = 'rev_phAveragesByPairings_ChAT_outcomeOnly';
ensureFigure(saveName, 1);
window = [2 5];
figsize = [2.4 2.3];
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & periReversalTrials, csPlusTrials & rewardTrials & extraReversalTrials}, 1,...
'FluorDataField', 'ZS', 'window', window, 'linespec', {'b', 'k'}); %high value, reward
% legend(hl, {['\bf\color{cyan}'], '\bf\color{magenta}'}, 'Location', 'northoutside', 'FontSize', 16, 'Box', 'off'); 
ylabel('\bf\fontsize{20}\color{green} Cholinergic'); set(gca, 'XLim', window, 'YLim', [-.5 3]);
xlabel('Time from reward (s)');
addStimulusPatch(gca, [2.9 3.1], '', patchHue);
set(gca, 'XTickLabel', {'-1', '0', '1', '2'});
formatFigureTalk(figsize);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.meta']));   
end
    
% just crop legend in DAT version in powerpoint
saveName = 'rev_phAveragesByPairings_DAT_outcomeOnly';
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & periReversalTrials, csPlusTrials & rewardTrials & extraReversalTrials}, 2,...
'FluorDataField', 'ZS', 'window', window, 'linespec', {'b', 'k'}); %high value, reward
% legend(hl, {['\bf\color{cyan} Initial trials'], '\bf\color{magenta} Others'}, 'Location', 'northoutside', 'FontSize', 16, 'Box', 'off'); 
set(gca, 'XLim', window, 'YLim', [-.5 3]);
ylabel('\bf\fontsize{20}\color{red} Dopaminergic');
xlabel('Time from reward (s)');
addStimulusPatch(gca, [2.9 3.1], '', patchHue);
set(gca, 'XTickLabel', {'-1', '0', '1', '2'});
formatFigureTalk(figsize);

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.meta']));   
end






















%% lets find optimal windows for CS for ACh and Dopamine

csPlus_special = csPlusTrials & excludeFirstTrials & TE.csLicksAll.rate > threshRate;
csMinus_special = csMinusTrials & excludeFirstTrials & TE.csLicksAll.rate <= threshRate;

ensureFigure('find_cue_windows', 1);
subplot(2,1,1);
[ha, hl] = phPlotAverageFromTE(TE, {csPlus_special, csMinus_special}, 1,...
'FluorDataField', 'ZS', 'window', [-.5, 1.5], 'linespec', {'g', 'r'}); %high value, reward
set(gca, 'XMinorTick', 'on', 'XGrid', 'on', 'XMinorGrid', 'on', 'GridColor', [0.3 0.3 0.3], 'MinorGridColor', [0.3 0.3 0.3], 'TickDir', 'both');
title('ACh');

subplot(2,1,2);
[ha, hl] = phPlotAverageFromTE(TE, {csPlus_special, csMinus_special}, 2,...
'FluorDataField', 'ZS', 'window', [-.5, 1.5], 'linespec', {'g', 'r'}); %high value, reward
set(gca, 'XMinorTick', 'on', 'XGrid', 'on', 'XMinorGrid', 'on', 'GridColor', [0.3 0.3 0.3], 'MinorGridColor', [0.3 0.3 0.3], 'TickDir', 'both');
title('Dopamine');

% set window to be 0.25 - 1 second after cue in LNL_odor_v2_pav_rev_AS
    


    







    
    
    
    