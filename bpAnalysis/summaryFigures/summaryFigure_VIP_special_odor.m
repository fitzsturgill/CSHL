function session = summaryFigure_VIP_special_odor(session, fCh)
    % one off analysis script for Adam's grant 
    
    
    zeroField = 'DeliverStimulus';

    lickHistSpecs = [-3, 6, 0.25]; % start, stop and bin size of lick histograms
    figName = 'VIP_special_odor';
    fig = ensureFigure(figName, 1); % ensure figure, erase if preexisting
    fig=mcPortraitFigSetup(fig);
    
    
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
    columns = 3;
    nAxes = 9;
    
    % layout axes in grid
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);
    
    
%%  Lick Rasters, etc.   
    bpLickRaster(SessionData, 1, 1, zeroField, '', hAx(1)); title('Cued: Reward');
    bpLickRaster(SessionData, 2, 2, zeroField, '', hAx(2)); title('Cued: Omit');
    bpLickRaster(SessionData, 3, 3, zeroField, '', hAx(3)); title('Uncued: Reward');
    bpLickHist(SessionData, [1 3], [1 3], lickHistSpecs,...
        'DeliverStimulus', 'PreCsRecording', 'PostUsRecording', {'m', 'r'}, [], hAx(7));
    title('Lick Hist');
    set(hAx([1:3 7]), 'XLim', [lickHistSpecs(1), lickHistSpecs(2)]); 
    
%% Photometry Rasters
    phSessionRaster(session.analysis, 1, 1, fCh, 'ax', hAx(4), 'lookupFactor', 4); ylabel('Trials');
    phSessionRaster(session.analysis, 2, 2, fCh, 'ax', hAx(5), 'lookupFactor', 4); 
    phSessionRaster(session.analysis, 3, 3, fCh, 'ax', hAx(6), 'lookupFactor', 4);      
       

%% Averages

    [~, hl] = phPlotSessionAverage(session.analysis, [1 2 3], [1 2 3], fCh, 'ax', hAx(8), 'linespec', {'m', 'k', 'r'});
    legend(hl, {'cued + R', 'omit', 'uncued R'}, 'Location', 'northwest'); legend('boxoff');
    xlabel('time (s)'); ylabel('dF/F');    
    
    set(hAx([4:6 8]), 'XLim', [-3, 6]);
    
%     set(hAx(10:12), 'YLim', [0.97 1.03]);
    
    
%% Saving
    saveas(fig, [fig_title '.fig']); %save as matlab fig
%     saveas(fig, [fig_title '.pdf']); %save as pdf
    disp('*** SummaryFigure complete ***');
    
    
    