% CuedOutcome_QPE_examples

%% Snippet to make lick histogram for graded value task from DC_35: (or DC_26?)

 % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-6 4], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
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