function session = SO_varyRewardSize_delay_summaryFigure(session, fCh)
% make and save summary figure for the SOv2_Training_FS protocol
% 2/22/16
% REPLACES TRIALOUTCOME VALUES WITH USOUTCOME VALUES!!!!!! TO PROPERLY USE
% bpFilterTrials

% assumes that 1) you've loaded session using session = bpLoadSession 

    if nargin < 2
        fCh = 1;
    end
    zeroField = 'DeliverStimulus';

    lickHistSpecs = [-3, 4, 0.25]; % start, stop and bin size of lick histograms
    figName = 'SO_varyRewardSize_delay_summaryFigure';
    fig = ensureFigure(figName, 1); % ensure figure, erase if preexisting
    fig=mcLandscapeFigSetup(fig);
    
    
    if ~isfield(session.SessionData, 'demod')
        session.SessionData = demodulateSession(session.SessionData);        
    end
    
    if ~isfield(session, 'analysis') || ~isfield(session.analysis, 'Photometry')
        session.analysis = processAnalysis_Photometry(session.SessionData, [], 'zeroField', zeroField);
    end

    %%
    SessionData = session.SessionData;


    
%% Define the axes matrix positions on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    matpos_title = [0 0 1 .1];
    matpos_big = [0 .1 1 .9];
%     matpos_Ph = [0 .5 1 .5];
    params.cellmargin = [.05 .05 0.05 0.05];    
    
%% Figure Title
    [~, fig_title, ~] = fileparts(session.filename); % session name
    title_ax = textAxes(fig, fig_title, 12);
    params.matpos = matpos_title;
    setaxesOnaxesmatrix(title_ax, 1, 1, 1, params, fig);
    
%% Setup Axes Matrix
%     trialTypes = [1 2; 1 3; 1 4; 2 2; 2 3; 2 4];
    params.matpos = matpos_big;
    rows = 3;
    columns = 4;
    nAxes = 12;
    
    % layout axes in grid
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);
    
    
%% first row:  Lick Rasters    
    bpLickRaster(SessionData, 1, 1, zeroField, '', hAx(1)); title('Cued: Small');
    bpLickRaster(SessionData, 2, 2, zeroField, '', hAx(2)); title('Cued: Big');
    bpLickRaster(SessionData, 4, 1, zeroField, '', hAx(3)); title('Uncued: Small');
    bpLickRaster(SessionData, 5, 2, zeroField, '', hAx(4)); title('Uncued: Big');    
    set(hAx(1:4), 'XLim', [lickHistSpecs(1), lickHistSpecs(2)]); 
    
%% second row: Photometry Rasters ( plus 1 more panel)   
    phSessionRaster(session.analysis, 1, 1, fCh, 'ax', hAx(5)); ylabel('Trials');
    phSessionRaster(session.analysis, 2, 2, fCh, 'ax', hAx(6)); 
    phSessionRaster(session.analysis, 4, 1, fCh, 'ax', hAx(7));
    phSessionRaster(session.analysis, 5, 2, fCh, 'ax', hAx(8));     
    phSessionRaster(session.analysis, 3, 3, fCh, 'ax', hAx(9)); title('Cued: Omit');         
       

%% third row: Photometry Averages

    [~, hl] = phPlotSessionAverage(session.analysis, [1 2 3], [1 2 3], fCh, 'ax', hAx(10), 'linespec', {'g', 'k', 'r'});
    legend(hl, {'small', 'big', 'omit'}, 'Location', 'northwest'); legend('boxoff');
    title('Cued');
    
    [~, hl] = phPlotSessionAverage(session.analysis, [4 5 6], [1 2 3], fCh, 'ax', hAx(11), 'linespec', {'g', 'k', 'r'});
    legend(hl, {'small', 'big', 'omit (valve click only)'}, 'Location', 'northwest'); legend('boxoff');
    title('Uncued');
    
    [~, hl] = phPlotSessionAverage(session.analysis, [1 4 2 5], [1 1 2 2], fCh, 'ax', hAx(12), 'linespec', {'c', 'b', 'm', 'r'});    
    legend(hl, {'cued Sm', 'uncued Sm', 'cued Bg', 'uncued Bg'}, 'Location', 'northwest'); legend('boxoff');
    xlabel('time (s)'); title('Cued vs Uncued');
    
%     set(hAx(10:12), 'YLim', [0.97 1.03]);
    
    
%% Saving
    saveas(fig, [fig_title '.fig']); %save as matlab fig
%     saveas(fig, [fig_title '.pdf']); %save as pdf
    disp('*** SummaryFigure complete ***');
    
    
    