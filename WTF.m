%% averages

% -3 6
window = [-1 5];
PhotometryField = 'Photometry';
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);


    saveName = ['Averages_405test'];
    fig = ensureFigure(saveName, 1);
    mcPortraitFigSetup(fig);

axes; hold on;    
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1]), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'g'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');

    

        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1]), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'m'}, 'PhotometryField', 'Photometry_ch2'); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        set(gca, 'XLim', window); ylabel('Fluor ZS');
