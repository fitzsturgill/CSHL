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
    ylabel('Fluorescence (Z Score)'); xlabel('Time from cue (s)');
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
%% Lick and Ph rasters from DC_20
CLimFactor = 2;
CLimFactor2 = 2.5;
trialRange = [170 350]; % global trial range
reversals = find(diff(TE.BlockNumber(csPlusTrials)));

    saveName = ['researchStatement_reversals_phRasters_dualChannel'];
    h=ensureFigure(saveName, 1);
%     mcPortraitFigSetup(h);
    

%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,3,1); 
    eventRasterFromTE(TE, csPlusTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
%     title('CS+'); ylabel('trial number');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', trialRange, 'TickDir', 'Out');
    set(gca, 'FontSize', 10)
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    set(gca, 'XLim', [-2 6]);
    
    subplot(1,3,2); phRasterFromTE(TE, csPlusTrials, 1, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 10, 'TickDir', 'Out')
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines        
            set(gca, 'YLim', trialRange);set(gca, 'XLim', [-2 6]);
     
    subplot(1,3,3); phRasterFromTE(TE, csPlusTrials, 2, 'CLimFactor', CLimFactor2, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 10, 'TickDir', 'Out')
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines       
        set(gca, 'YLim', trialRange);set(gca, 'XLim', [-2 6]);
  


    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
    end

    
    
    
    
    