% SO_RewardPunish_figure

DB = dbLoadExperiment('SO_RewardPunish_odor');
saveOn = 1;
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);

%% averages from ChAT_22

animal = 'ChAT_22';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups

%% photometry averages, zscored
%     ylim = [-2 8];
    window = [-5 4];
    fdField = 'ZS';
    saveName = 'ChAT_22_phAverages_appetitive';
    ensureFigure(saveName, 1);
    
    
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5]), 1,...
        'PhotometryField', 'Photometry', 'FluorDataField', fdField, 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'b','k','c'});     
    set(gca, 'XLim', window);
    ylabel('(\fontsize{12}\sigma\fontsize{8}-baseline)'); xlabel('Time from reinforcement (s)');
    addStimulusPatch(gca, [-2 -1], '', [0.8 0.8 0.8]); addStimulusPatch(gca, [-0.1 0.1], '', [0.8 0.8 0.8]);
%     legend(hl, {'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    formatFigurePublish('size', [1.5 1]);    
    
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end    