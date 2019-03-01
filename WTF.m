%% averages

% Define the positions of different axes matrices on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    params = struct();
    matpos_title = [0 0 1 .1];
    matpos_avgs = [0 .1 1 0.85];    
    params.cellmargin = [.05 .05 0.05 0.05];    
    params.figmargin = [0.05 0.05 0.025 0.025];


    
    
    % -3 6
window = [-1 5];
PhotometryField = 'Photometry';
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);

% rewarding subset
showTheseR = find(Odor2Valve1Trials & rewardTrials);
[~, rixR] = sort(rand(size(showTheseR)));
% aversive subset
showTheseA = find(Odor2Valve2Trials & punishTrials);
[~, rixA] = sort(rand(size(showTheseA)));


for channel = channels
    saveName = ['Averages_ch' num2str(channel)];
    fig = ensureFigure(saveName, 1);
    mcPortraitFigSetup(fig);

    params.matpos = matpos_title;
    [ax, ~] = textAxes(fig, sprintf('%s , Channel %d', subjectName, channel), 14);
    setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);    
    
%     try
        ax = subplot(3,2,1);
        plot(TE.Photometry.xData, TE.Photometry.data(channel).ZS(showTheseR(rixR(1:10)), :)', 'LineWidth', 0.25); set(gca, 'XLim', window);
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('Fluor. (ZS)'); title('rewarding');
        disp('hello');
%     catch
%     end
    
    try
        ax(2) = subplot(3,2,2);
        plot(TE.Photometry.xData, TE.Photometry.data(channel).ZS(showTheseA(rixA(1:10)), :)', 'LineWidth', 0.25); set(gca, 'XLim', window);
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('Fluor. (ZS)'); title('aversive');
    catch
    end    

    ax(3) = subplot(3,2,3);
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5]), channel, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k', 'c'}, 'PhotometryField', PhotometryField); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
    set(gca, 'XLim', window); ylabel('Fluor ZS');
    
    try
        ax(4) = subplot(3,2,4);
        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([3 4 6]), channel, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k', 'm'}, 'PhotometryField', PhotometryField); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        set(gca, 'XLim', window); ylabel('Fluor ZS');
    catch
    end
    
    try
        ax(5) = subplot(3,2,5);
        varargin = {'window', window, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k','c'}};
        axh = [];
        [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 2 5]), 'Port1In', varargin{:});  
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('licks (1/s)'); xlabel('time from cue(s)');  set(gca, 'XLim', window);
        legend(hl, {'cued reward', 'cued ommission', 'uncued reward'}, 'Box', 'off', 'Location', 'best'); 
    catch
    end
    
    try
        ax(6) = subplot(3,2,6);
        varargin = {'window', window, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k','m'}};
        axh = [];
        [ha, hl] = plotEventAverageFromTE(TE, trialsByType([3 4 6]), 'Port1In', varargin{:});  
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('licks (1/s)'); xlabel('time from cue(s)');  set(gca, 'XLim', window); legend(hl, {'cued punish', 'cued ommission', 'uncued punish'}, 'Box', 'off', 'Location', 'best');      
    catch
    end
    

    
    params.matpos = matpos_avgs;
    setaxesOnaxesmatrix(ax, 3, 2, 1:6, params, fig);     
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end    
end