% CuedOutcome_QPE_examples
DB = dbLoadExperiment('cuedOutcome');
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);
saveOn = 1;


% photometryField = 'Photometry';
% fdField = 'ZS';


%% Lick histogram for graded value task
animal = 'ChAT_39';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups
 % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-6 4], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};

    ensureFigure('Lick_Hist_CCN', 1); axes;

    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 4 7]), 'Port1In', varargin{:});
    set(gca, 'XLim', [-6 4], 'YLim', [-2 15]);
    addStimulusPatch(gca, [-3 -2 10], '', [0.7 0.7 0.7]) 
    addStimulusPatch(gca, [-0.1 0.1], '');
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'FontSize', 6, 'Interpreter', 'tex', 'Position', [0.25 0.6 0.3 0.4]); legend('boxoff');    
    ylabel('Licks (Hz)'); xlabel('Time from reinforcement (s)');
    
    formatFigurePublish('size', [2 1.1]);
    if saveOn 
        export_fig(fullfile(savepath, saveName), '-eps');
    end    
    
