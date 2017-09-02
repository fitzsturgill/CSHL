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
    axh = [];
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
        saveas(gcf, fullfile(savepath, 'Ph_Hist_Dopamine_CCN.epsc'));           
    end
    
    
%% Snippet to make graded value averages from ChAT_42

    ensureFigure('GradedValue_Cue', 1); 
    axes(); [ha, hla] = phPlotAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 1,...
        'window', [-4 0], 'linespec', {'b', 'r', 'k'}, 'FluorDataField', 'ZS');
    set(gca, 'XLim', [-4 0]);
    formatFigureGRC;
    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'GradedValue_Cue.fig'));
        saveas(gcf, fullfile(savepath, 'GradedValue_Cue.jpg'));    
        saveas(gcf, fullfile(savepath, 'GradedValue_Cue.epsc'));           
    end
    
%% Snippet to make graded value lick averages from ChAT_42

 % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-6 0], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};
    axh = [];
    ensureFigure('GradedValue_Cue_Lick', 1); axes('FontSize', 12, 'LineWidth', 1);
    plotEventAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 'Port1In', varargin{:});
    set(gca, 'XLim', [-6 0]);
    formatFigureGRC;    
    ylabel('Licks (Hz)'); xlabel('Time from cue (s)');
    
    if saveOn    
        saveas(gcf, fullfile(savepath, 'GradedValue_Cue_Lick.fig'));
        saveas(gcf, fullfile(savepath, 'GradedValue_Cue_Lick.jpg'));    
        saveas(gcf, fullfile(savepath, 'GradedValue_Cue_Lick.epsc'));           
    end    
    
    
    
    
    
    