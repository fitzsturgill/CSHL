% CCN Talk Script

%% desktop
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\Adam\';
saveOn = 1;

%% laptop
savepath = 'C:\Users\Fitz\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\Adam\';
saveOn = 1;
%% Snippet to make lick histogram for graded value task from DC_35: (or DC_26?)

 % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-6 4], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...a
        'linespec', {'b', 'r', 'g'}};

    ensureFigure('Lick_Hist_CCN', 1); axes;

    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 4 7]), 'Port1In', varargin{:});
    set(hl, 'LineWidth', 2);
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    set(gca, 'XLim', [-6 4], 'YLim', [-2 15]);
    addStimulusPatch(gca, [-3 -2 10], '', [0.7 0.7 0.7]) 
    addStimulusPatch(gca, [-0.1 0.1], '');
    ylabel('Licks (Hz)'); xlabel('Time from reinforcement (s)');
    formatFigureTalk([4 3]);
    if saveOn    
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN.fig'));
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN.jpg'));    
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN.emf'));
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN'), 'epsc');        
    end
    
%% Us and Cs scatter plots
cuedOutcome_pooledAnalysis_script2;
  
%% Snippet to make dopamine photometry averages from DC_26:

    ensureFigure('Ph_Hist_Dopamine_CCN', 1); axes('FontSize', 12, 'LineWidth', 1); 
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1, 'FluorDataField', 'ZS', 'linespec', {'b', 'r', 'g'}); % reward, varying degrees of expectation
    set(hl, 'LineWidth', 2);
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');    
    set(gca, 'XLim', [-5 4], 'YLim', [-2 16]); ylabel('Z Score'); xlabel('time from reinforcement (s)');
    addStimulusPatch(gca, [-3 -2 10], '', [0.7 0.7 0.7]) 
    addStimulusPatch(gca, [-0.1 0.1], '');
    ylabel('Dopamine (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from reinforcement (s)');
    formatFigureTalk([4 3]);
    if saveOn    
        saveas(gcf, fullfile(savepath, 'Ph_Hist_Dopamine_CCN.fig'));
        saveas(gcf, fullfile(savepath, 'Ph_Hist_Dopamine_CCN.jpg'));    
        saveas(gcf, fullfile(savepath, 'Ph_Hist_Dopamine_CCN.emf'));           
        saveas(gcf, fullfile(savepath, 'Ph_Hist_Dopamine_CCN.pdf'));
    end
    
    
%% Snippet to make cuedOutcome cue response photometry histogram from ChAT_42

    ensureFigure('CuedOutcome_Cue', 1); 
    [ha, hl] = phPlotAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 1,...
        'window', [-4 3], 'linespec', {'b', 'r', 'g'}, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue);
    set(hl, 'LineWidth', 2);
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    set(gca, 'XLim', [-4 3]);
    addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7]) 
    ylabel('Cholinergic (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from cue (s)');
    formatFigureTalk([4 3]);    

    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue.fig'));
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue.jpg'));    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue.emf'));           
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue.svg'));            
    end
    
%% Snippet to make cuedOutcome cue response Lick histogram from ChAT_42
 % cue types
    varargin = {'window', [-4 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};
    ensureFigure('CuedOutcome_Cue_Lick', 1);
    [ha, hl] = plotEventAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 'Port1In', varargin{:});
    set(hl, 'LineWidth', 2);
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    set(gca, 'XLim', [-4 3], 'YLim', [0 5]);
    addStimulusPatch(gca, [0 1], 'odor', [0.7 0.7 0.7]) 
    ylabel('Licks (Hz)'); xlabel('Time from cue (s)');
    formatFigureTalk([4 2]);        
    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_Lick.fig'));
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_Lick.jpg'));    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_Lick.emf'));           
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_Lick.svg'));           
    end
    
%% Snippet to make combined cue PSTH from ChAT_42 (licking and photometry)

    ensureFigure('CuedOutcome_Cue_combined', 1); 
        varargin = {'window', [-4 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};
    lickAvg = eventAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 'Port1In', varargin{:});
    ax = axes; hold on; yyaxis right
    ll = plot(lickAvg.xData, lickAvg.Avg(1,:), '--k', 'LineWidth', 2);
    plot(lickAvg.xData, lickAvg.Avg(1,:), '--b', lickAvg.xData, lickAvg.Avg(2,:), '--r', lickAvg.xData, lickAvg.Avg(3,:), '--g', 'LineWidth', 2);
    ax.YColor = [0 0 0]; ylabel('Lick rate (Hz)');
    yyaxis left; ax.YColor = [0 0 0];
    [ha, hl] = phPlotAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 1,...
        'window', [-4 3], 'linespec', {'b', 'r', 'g'}, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue);
    set(hl, 'LineWidth', 2);    
    set(gca, 'XLim', [-4 3], 'YLim', [-0.2 1.5]);
    addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7]) 
    legend([hl ll(1)], {'\color{blue} high value', '\color{red} low value', '\color{green} uncued', ...
        'licking'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    ylabel('Cholinergic (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from cue (s)');
    formatFigureTalk([4 3]);    

    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_combined.fig'));
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_combined.jpg'));    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Cue_combined.emf'));                
    end    
    
    
    %% Snippet to make outcome PSTH (reward condition from ChAT_42
    
    ensureFigure('CuedOutcome_Reward', 1); 
    [ha, hl] = phPlotAverageFromTE(TE, {trialsByType{1}, trialsByType{4}, trialsByType{7}}, 1,...
        'window', [-2 2], 'linespec', {'b', 'r', 'g'}, 'FluorDataField', 'ZS');
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
    set(hl, 'LineWidth', 2);
    set(gca, 'XLim', [-2 2], 'YLim', [-1.2 2.5]);
%     addStimulusPatch(gca, [-3 -2 10], '', [0.7 0.7 0.7]) 
    addStimulusPatch(gca, [-0.1 0.1], '');
    ylabel('Cholinergic (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from reinforcement (s)');
    formatFigureTalk([4 3]);    

    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward.fig'));
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward.jpg'));    
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward.emf'));           
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward.svg'));           
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
        saveas(gcf, fullfile(savepath, 'CuedOutcome_Reward_Lick.emf'));           
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
    set(gca, 'XLim', [-2 6], 'XTick', [0 3 6]); ylabel('Trial #');
    addStimulusBar(gca, [0 1 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [2.9 3.1 0], '', [0.3 0.3 0.3], 5);
    
    subplot(1,3,2); phRasterFromTE(TE, csPlusTrials & includedTrials & rewardTrials, 1, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-2 6]);
        set(gca, 'FontSize', 10, 'TickDir', 'Out')
    line(repmat([-2; 6], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % reversal lines        
        set(gca, 'XLim', [-2 6], 'YTick', [], 'XTick', [0 3 6], 'Clipping', 'off'); xlabel('Time from cue (s)');
    addStimulusBar(gca, [0 1 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [2.9 3.1 0], '', [0.3 0.3 0.3], 5); 
     
    subplot(1,3,3); phRasterFromTE(TE, csPlusTrials & includedTrials & rewardTrials, 2, 'CLimFactor', CLimFactor2, 'trialNumbering', 'consecutive', 'window', [-2 6]);
        set(gca, 'FontSize', 10, 'TickDir', 'Out')
    line(repmat([-2; 6], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % reversal lines       
        set(gca, 'XLim', [-2 6], 'Clipping', 'off');
    set(gca, 'XLim', [-2 6], 'YTick', [], 'XTick', [0 3 6]);
    addStimulusBar(gca, [0 1 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [2.9 3.1 0], '', [0.3 0.3 0.3], 5);    
    formatFigureTalk([4.5 3]);        
  


    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.emf']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));           
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
    formatFigurePoster([6 4], '', 20);

    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.emf']));   
    end
    
%% sorting by behavior (doesn't really work well)


excludeFirstNTrials = 20;
excludeFirstTrials = TE.trialNumber > excludeFirstNTrials;
threshPercentile = 0.15;

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
trialMatrix2 = trialMatrix;
for counter = 1:nSessions
    trialsThisSession = TE.sessionIndex == counter;
    thisReverse = find(TE.BlockChange & trialsThisSession);
    postReversalTrials = trialsThisSession; postReversalTrials(1:thisReverse - 1) = 0;
    pairingsThisReverse = postReversalTrials & csPlusTrials & rewardTrials;    
    trialSearch = find(pairingsThisReverse);
    threshTrial = trialSearch(rewardPairingsThresh);
    firstTrials = postReversalTrials; firstTrials(threshTrial:end) = 0;
    trialMatrix(:,counter) = firstTrials;
    threshTrial2 = trialSearch(rewardPairingsThresh * 2 + 1);
    nextTrials = postReversalTrials; nextTrials([1:threshTrial - 1 threshTrial2:end]) = 0;
    trialMatrix2(:,counter) = nextTrials;
end
periReversalTrials = any(trialMatrix, 2);
extraReversalTrials = ~periReversalTrials;
% extraReversalTrials = any(trialMatrix2, 2);
% find the next 20

% extraReversalTrials = extraReversalTrials(

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
    saveas(gcf, fullfile(savepath, saveName), 'epsc');
    saveas(gcf, fullfile(savepath, saveName), 'fig');   
    saveas(gcf, fullfile(savepath, saveName), 'meta');   
    saveas(gcf, fullfile(savepath, saveName), 'svg');      
    saveas(gcf, fullfile(savepath, saveName), 'pdf');          
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
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
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
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
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
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
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




%% Snippet to make uncued reward and punishment responses from ChAT_26
    saveName = ['Reward_ChAT'];
    ensureFigure(saveName, 1); 
    [ha, hl] = phPlotAverageFromTE(TE, rewardTrials & uncuedTrials, 1,...
        'window', [-1 3], 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'cmap', [171 55 214]/256);
    set(hl, 'LineWidth', 2);
    set(gca, 'XLim', [-1 3], 'YLim', [-1 2], 'YTick', [-1 0 1], 'XTick', []);
    addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7]) 
    ylabel('Fluor. (\fontsize{20}\sigma\fontsize{16}-baseline)'); %xlabel('Time from reinforcement (s)');
    formatFigureTalk([3.5 2.5]);    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
end    
%% Snippet to make uncued reward and punishment responses from ChAT_26
    saveName = ['Punish_ChAT'];
    ensureFigure(saveName, 1); 
    [ha, hl] = phPlotAverageFromTE(TE, punishTrials & uncuedTrials, 1,...
        'window', [-1 3], 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'cmap', [171 55 214]/256);
    set(hl, 'LineWidth', 2);
    set(gca, 'XLim', [-1 3], 'YLim', [-1 2], 'YTick', [-1 0 1]);
    addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7]) 
    ylabel('Fluor. (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from reinforcement (s)');
    formatFigureTalk([3.5 2.2]);    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
end    
    
%% Snippet to make uncued reward and punishment responses from DC_26
    saveName = ['Reward_DAT'];
    ensureFigure(saveName, 1); 
    [ha, hl] = phPlotAverageFromTE(TE, rewardTrials & uncuedTrials, 1,...
        'window', [-1 3], 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'cmap', [237 125 49]/256);
    set(hl, 'LineWidth', 2);
    set(gca, 'XLim', [-1 3], 'YLim', [-5 16], 'XTick', []);
    addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7]) 
    ylabel('Fluor. (\fontsize{20}\sigma\fontsize{16}-baseline)');% xlabel('Time from reinforcement (s)');
    formatFigureTalk([3.5 2.2]);    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
end    
%% Snippet to make uncued reward and punishment responses from DC_26
    saveName = ['Punish_DAT'];
    ensureFigure(saveName, 1); 
    [ha, hl] = phPlotAverageFromTE(TE, punishTrials & uncuedTrials, 1,...
        'window', [-1 3], 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'cmap', [237 125 49]/256);
    set(hl, 'LineWidth', 2);
    set(gca, 'XLim', [-1 3], 'YLim', [-3 3]);%, 'YTick', [-1 0 1]);
    addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7]) 
    ylabel('Fluor. (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from reinforcement (s)');
    formatFigureTalk([3.5 2.2]);    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
end    


%% snippet from cuedOutcome_pooledAnalysis_grandAverages
% plot grand average photometry responses for reward condition, don't normalize (not clear how you would)
% 3rd dimension - first chunk = high value, second = low value, third =
% uncued

saveOn = 1;
saveName = 'Full_grandAverage_SEM_rewardedOutcomes';
ensureFigure(saveName, 1); %axes; set(gca, 'YLim', [-0.2 0.8]);
% addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7]);

cmap = [0 0 1; 1 0 0; 0 1 0];
nSessions = size(avgData.full.licks.avg, 1);
xData_licks_full = squeeze(avgData.full.licks.xData(1, :, 1));
xData_ph_full = squeeze(avgData.full.photometry.xData(1, :, 1));
fullLicks_avg = squeeze(mean(avgData.full.licks.avg));
fullPh_avg = squeeze(mean(avgData.full.photometry.avg));
fullLicks_sem = squeeze(std(avgData.full.licks.avg)) / sqrt(nSessions);
fullPh_sem = squeeze(std(avgData.full.photometry.avg)) / sqrt(nSessions);

% first axes column- reward, second axes column- punishment, 3rd, - neutral
axes; set(gca, 'YLim', [-0.25 2.25]); addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7]); addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7]);
[hl, hp] = boundedline(xData_ph_full', fullPh_avg(:,[1 4 7]), permute(fullPh_sem(:,[1 4 7]), [1 3 2]), 'cmap', cmap);
set(hl, 'LineWidth', 2);
legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 16, 'Interpreter', 'tex'); legend('boxoff');
set(gca, 'XLim', [-5 3]);
ylabel('Cholinergic (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from reinforcement (s)');

formatFigureTalk([4 3]);
if saveOnh
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
    saveas(gcf, fullfile(savepath, [saveName]), 'epsc');       
end    
    
    