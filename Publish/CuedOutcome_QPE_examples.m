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
    
    
%% combined licking and photometry averages from ChAT_42


    load('Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_42\TE.mat');
    cuedOutcome_Conditions;
    saveName = 'CuedOutcome_example_averages_combined_ChAT_42';
    window = [-1.5 5];
    ensureFigure(saveName, 1); 
        varargin = {'window', [window(1) 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};
    lickAvg = eventAverageFromTE(TE, {highValueTrials & rewardTrials, lowValueTrials & rewardTrials, uncuedTrials & rewardTrials}, 'Port1In', varargin{:});
    ax = axes; hold on; yyaxis right
    ll = plot(lickAvg.xData, lickAvg.Avg(1,:), '--k');
    plot(lickAvg.xData, lickAvg.Avg(1,:), '--b', lickAvg.xData, lickAvg.Avg(2,:), '--r', lickAvg.xData, lickAvg.Avg(3,:), '--g');
    ax.YColor = [0 0 0]; ylabel('Lick rate (Hz)');
    set(gca, 'YLim', [-1 7]);
    yyaxis left; ax.YColor = [0 0 0]; 
    [ha, hl] = phPlotAverageFromTE(TE, {highValueTrials & rewardTrials, lowValueTrials & rewardTrials, uncuedTrials & rewardTrials}, 1,...
        'window', window, 'linespec', {'b', 'r', 'g'}, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue, 'alpha', 0);
%     set(hl, 'LineWidth', 2);    
    set(gca, 'XLim', window, 'YLim', [-1 3]);
    addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5) 
    addStimulusPatch(gca, [2.9 3.1], '', [0.8 0.8 0.8], 0.5) 
%     legend([hl ll(1)], {'\color{blue} high value', '\color{red} low value', '\color{green} uncued', ...
%         'licking'}, 'Location', 'northwest', 'FontSize', 20, 'Interpreter', 'tex'); legend('boxoff');
    ylabel('Fluor. (\sigma-baseline)'); xlabel('Time from cue (s)');
    formatFigurePublish('size', [1.7 1.5]);    

if saveOn
    export_fig(fullfile(savepath, saveName), '-eps');
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
%     saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    

    
