function session = summaryFigure_SO_rewardPunish_odor(session, fCh)


    zeroField = 'DeliverStimulus';
    try
        PunishOn = session.SessionData.TrialSettings(1).GUI.PunishOn;
    catch
        PunishOn = 1;
    end

%     lickHistSpecs = [-3, 5, 0.25]; % start, stop and bin size of lick histograms
%     lickHistSpecs = 0.25;
    figName = ['SO_varyRewardSize_delay_summaryFigure_ch' num2str(fCh)];
    fig = ensureFigure(figName, 1); % ensure figure, erase if preexisting
    fig=mcLandscapeFigSetup(fig);
    
    % pull out relevent info from trial settings
    ts = session.SessionData.TrialSettings(1);
    
    x1 = -ts.PreCsRecording + 1; % cut off first second of recording due to transient in photometry
    x2 = ts.GUI.OdorTime + ts.GUI.Delay + ts.PostUsRecording;
    odorTime =  [0 ts.GUI.OdorTime];
    UsTime = [(ts.GUI.OdorTime + ts.GUI.Delay - 0.05) (ts.GUI.OdorTime + ts.GUI.Delay + 0.05)];
    lickHistSpecs = [x1 x2 0.25];
    
    
    if ~isfield(session.SessionData, 'demod')
        session.SessionData = demodulateSession(session.SessionData);        
    end
    
    if ~isfield(session, 'analysis') || ~isfield(session.analysis, 'Photometry')
        session.analysis.Photometry = processAnalysis_Photometry(session,'zeroField', zeroField, 'downsample', 10);
    end

    %%
    SessionData = session.SessionData;


    
%% Define the axes matrix positions on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    matpos_title = [0 0 1 .1];
%     matpos_big = [0 .1 1 .9];
    matpos_lickRaster = [0 0.1 2/5 0.6 * 0.9]; % 2/3 * 0.9,  2/3 of fig height discounting height of title axis
    matpos_phRaster = [2/5 0.1 3/5 0.6 * 0.9];
    matpos_avgs = [0 (0.6 * 0.9 + .1) 1 0.4 * 0.9];
%     matpos_Ph = [0 .5 1 .5];
    params.cellmargin = [.05 .05 0.05 0.05];    
    
%% Figure Title
    [~, fig_title, ~] = fileparts(session.filename); % session name
    fig_title = [fig_title '_ch' num2str(fCh)];
    title_ax = textAxes(fig, fig_title, 12);
    params.matpos = matpos_title;
    setaxesOnaxesmatrix(title_ax, 1, 1, 1, params, fig);
    
%% Setup Lick Rasters Axes Matrix
%     trialTypes = [1 2; 1 3; 1 4; 2 2; 2 3; 2 4];
    params.matpos = matpos_lickRaster;
    rows = 2;
    columns = 2;
    nAxes = 4;
    
    % layout axes in grid
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);
    
    
%% Lick Rasters    
% outcomes  1- reward, 2 - punish, 3- omission
    bpLickRaster(SessionData, 1, 1, zeroField, '', hAx(1)); title('Cued'); ylabel('Reward Odor');
    stimBars(gca, 0.1, 0, odorTime, UsTime);
    bpLickRaster(SessionData, 5, 1, zeroField, '', hAx(2)); title('Uncued');stimBars(gca, 0.1, 0, odorTime, UsTime);
    if PunishOn
        bpLickRaster(SessionData, 3, 2, zeroField, '', hAx(3)); ylabel('Punish Odor'); stimBars(hAx(3), 0.1, 1, odorTime, UsTime);    
        bpLickRaster(SessionData, 6, 2, zeroField, '', hAx(4));  stimBars(hAx(4), 0.1, 1, odorTime, UsTime);
    end
    set(hAx(1:4), 'XLim', [x1, x2]); 
    
%% Setup Ph Rasters Axes Matrix
%     trialTypes = [1 2; 1 3; 1 4; 2 2; 2 3; 2 4];
    params.matpos = matpos_phRaster;
    rows = 2;
    columns = 3;
    nAxes = 6;
    lookupFactor = 4; %4;
    
    % layout axes in grid
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);
    %% second row: Photometry Rasters ( plus 1 more panel)   
    phSessionRaster(session.analysis, 1, 1, fCh, 'ax', hAx(1), 'lookupFactor', lookupFactor); title('Cued: Us'); ylabel('Rewarded Odor');
    phSessionRaster(session.analysis, 2, 3, fCh, 'ax', hAx(2), 'lookupFactor', lookupFactor); title('Cued: Omit'); 
    phSessionRaster(session.analysis, 5, 1, fCh, 'ax', hAx(3), 'lookupFactor', lookupFactor); title('Uncued: Us');
    if PunishOn
        phSessionRaster(session.analysis, 3, 2, fCh, 'ax', hAx(4), 'lookupFactor', lookupFactor);     ylabel('Punished Odor');
        phSessionRaster(session.analysis, 4, 3, fCh, 'ax', hAx(5), 'lookupFactor', lookupFactor); 
        phSessionRaster(session.analysis, 6, 2, fCh, 'ax', hAx(6), 'lookupFactor', lookupFactor);
    end
    
    for i = 1:length(hAx)
        axes(hAx(i));
        yl = get(gca, 'YLim');
        yPos = abs(diff(yl)) * 0.1 + min(yl); 
        addStimulusBar(gca, [odorTime yPos], '', [1 1 1], 5);
        addStimulusBar(gca, [UsTime yPos], '', [1 1 1], 5);
    end
       

%% Setup Averages Axes Matrix
%     trialTypes = [1 2; 1 3; 1 4; 2 2; 2 3; 2 4];
    params.matpos = matpos_avgs;
    rows = 1;
    columns = 3;
    nAxes = 3;
    
    % layout axes in grid
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);

    %% third row: Averages
    bpLickHist(SessionData, [1 5], [1 1], lickHistSpecs,...
        'DeliverStimulus', 'PreCsRecording', 'PostUsRecording', {'c', 'b'}, [], hAx(1));
    stimBars(hAx(1), 0.9, 1, odorTime, UsTime);
    title('Lick Hist: Rewarded Odor'); xlabel('time (s)'); ylabel('Lick rate (1/s)');

    [~, hl] = phPlotSessionAverage(session.analysis, [1 2 5], [1 3 1], fCh, 'ax', hAx(2), 'linespec', {'c', 'k', 'b'});
    legend(hl, {'cued', 'omit', 'uncued'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    title('deltaF/F: Rewarded Odor'); ylabel('dF/F'); xlabel('time (s) from odor R'); 
    stimBars(hAx(2), 0.95, 0, odorTime, UsTime);
    
    if PunishOn
        [~, hl] = phPlotSessionAverage(session.analysis, [3 4 6], [2 3 2], fCh, 'ax', hAx(3), 'linespec', {'m', 'k', 'r'});
        legend(hl, {'cued', 'omit', 'uncued'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        title('deltaF/F: Punished Odor'); ylabel('dF/F'); xlabel('time (s) from odor P');   
    %     set(hAx(2:3), 'YLim', [-0.005 0.01]);
        stimBars(hAx(3), 0.95, 0, odorTime, UsTime);
    end
    
    set(hAx, 'XLim', [x1 x2]);
    
%% Saving
    saveas(fig, [fig_title '.fig']); %save as matlab fig
%     saveas(fig, [fig_title '.pdf']); %save as pdf
    disp('*** SummaryFigure complete ***');
end

function stimBars(ax, yPos, labelOn, odorTime, UsTime)
    if nargin < 3
        labelOn = 0;
    end
    yl = get(ax, 'YLim');
    yPos = abs(diff(yl)) * yPos + min(yl);
    cs = ''; us = '';
    if labelOn
        cs = 'odor';
        us = 'US';
    end
    
    addStimulusBar(ax, [odorTime yPos], cs, [0 0 0], 4); 
    addStimulusBar(ax, [UsTime yPos], us, [1 0 0], 4);
end
