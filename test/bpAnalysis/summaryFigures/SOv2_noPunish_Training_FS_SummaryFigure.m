function session = SOv2_noPunish_Training_FS_SummaryFigure(session, fCh)
% make and save summary figure for the SOv2_Training_FS protocol
% 2/22/16
% REPLACES TRIALOUTCOME VALUES WITH USOUTCOME VALUES!!!!!! TO PROPERLY USE
% bpFilterTrials

% assumes that 1) you've loaded session using session = bpLoadSession 

    if nargin < 2
        fCh = 1;
    end
    zeroField = 'DeliverStimulus';
    startField = 'ITI';
%     if noFluor
%         endField = {'PostUS', 'last'};
%     else
        endField = 'PostTrialRecording';
%     end
    lickHistSpecs = [-2, 4, 0.25]; % start, stop and bin size of lick histograms
    fig = ensureFigure('SOv2_noPunish_Training_FS_SummaryFigure', 1); % ensure figure, erase if preexisting
    fig=mcLandscapeFigSetup(fig);
    %% replace TrialOutcome with UsOutome values
    session.SessionData.TrialOutcome = session.SessionData.UsOutcome(1, 1:session.SessionData.nTrials) + 1;
    
    
    if ~isfield(session.SessionData, 'demod')
        session.SessionData = demodulateSession(session.SessionData);        
    end
    
    if ~isfield(session, 'analysis') || ~isfield(session.analysis, 'Photometry')
        session.analysis = processAnalysis_Photometry(session.SessionData, [], 'zeroField', 'DeliverStimulus', 'startField', 'PreTrialRecording');
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
    bpLickRaster(SessionData, 1, 2, zeroField, '', hAx(1)); title('toneA: Reward');
    bpLickRaster(SessionData, 2, 2, zeroField, '', hAx(2)); title('toneB: Reward');
    bpLickRaster(SessionData, 3, 2, zeroField, '', hAx(3)); title('toneA: Omit');
    bpLickRaster(SessionData, 4, 1, zeroField, '', hAx(4)); title('toneB: Omit');    
    set(hAx(1:4), 'XLim', [lickHistSpecs(1), lickHistSpecs(2)]); 
    
%% second row: Photometry Rasters    
    phSessionRaster(session.analysis, 1, 2, fCh, 'ax', hAx(5)); ylabel('Trials');
    phSessionRaster(session.analysis, 2, 2, fCh, 'ax', hAx(6)); 
    phSessionRaster(session.analysis, 3, 2, fCh, 'ax', hAx(7));
    phSessionRaster(session.analysis, 4, 1, fCh, 'ax', hAx(8));     
       

%% third row: Photometry Averages

    [~, hl] = phPlotSessionAverage(session.analysis, [1 1], [2 1], fCh, 'ax', hAx(9), 'linespec', {'g', 'k'});
    legend(hl, {'reward', 'omit'}, 'Location', 'northwest');
    title('tone A');
    
    [~, hl] = phPlotSessionAverage(session.analysis, [2 2], [2 1], fCh, 'ax', hAx(10), 'linespec', {'g', 'k'});
    legend(hl, {'reward', 'omit'}, 'Location', 'northwest');
    title('tone B');
    
    [~, hl] = phPlotSessionAverage(session.analysis, [1 2 3], [2 2 2], fCh, 'ax', hAx(11), 'linespec', {'c', 'm', 'b'});    
    legend(hl, {'toneA', 'toneB', 'uncued'}, 'Location', 'northwest');
    xlabel('time (s)'); title('reward');
    
    [~, hl] = phPlotSessionAverage(session.analysis, [1 2 4], [1 1 1], fCh, 'ax', hAx(12), 'linespec', {'c', 'm', 'b'});        
    legend(hl, {'toneA', 'toneB', 'uncued'}, 'Location', 'northwest');
    xlabel('time (s)'); title('omit');    
    
    
    
    
%% Saving
    saveas(fig, [fig_title '.fig']); %save as matlab fig
%     saveas(fig, [fig_title '.pdf']); %save as pdf
    disp('*** SummaryFigure complete ***');
    
    
    