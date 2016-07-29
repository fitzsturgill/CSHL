function session = SO_varyRewardSize_summaryFigure(session, fCh)

    if nargin < 2
        fCh = 1;
    end

    zeroField = 'DeliverStimulus';
    startField = 'PreCsRecording';
    lickHistSpecs = [-2, 4, 0.25]; % start, stop and bin size of lick histograms    
    fig = ensureFigure('SO_varyRewardSize_summaryFigure', 1); % ensure figure, erase if preexisting
    fig=mcPortraitFigSetup(fig);
    
    if ~isfield(session.SessionData, 'demod')
        session.SessionData = demodulateSession(session.SessionData);        
    end
    
    if ~isfield(session, 'analysis') || ~isfield(session.analysis, 'Photometry')
        session.analysis = processAnalysis_Photometry(session.SessionData, [], 'zeroField', 'DeliverStimulus');
    end
    
    SessionData = session.SessionData;
    
%% Define the axes matrix positions on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    matpos_title = [0 0 1 .1];
    matpos_Licks = [0 .1 1 .4];
    matpos_Ph = [0 .5 1 .5];
    params.cellmargin = [.05 .05 0.05 0.05];    
    
%% Figure Title
    [~, fig_title, ~] = fileparts(session.filename); % session name
    title_ax = textAxes(fig, fig_title, 12);
    params.matpos = matpos_title;
    setaxesOnaxesmatrix(title_ax, 1, 1, 1, params, fig);
    
%% Lick Rasters
%     trialTypes = [1 2; 1 3; 1 4; 2 2; 2 3; 2 4];
    params.matpos = matpos_Licks;
    rows = 2;
    columns = 4;
    nAxes = 8;
    % layout axes in grid
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);
    bpLickRaster(SessionData, 1, 1, zeroField, '', hAx(1)); ylabel('Short Tone'); title('Small Reward');
    bpLickRaster(SessionData, 2, 2, zeroField, '', hAx(2)); title('Big');
    bpLickRaster(SessionData, 3, 3, zeroField, '', hAx(3)); title('Omit');
    bpLickRaster(SessionData, 7, 1, zeroField, '', hAx(4)); title('Uncued Small');    
    bpLickRaster(SessionData, 4, 1, zeroField, '', hAx(5)); ylabel('Long Tone'); xlabel('time (s) from tone');
    bpLickRaster(SessionData, 5, 2, zeroField, '', hAx(6)); 
    bpLickRaster(SessionData, 6, 3, zeroField, '', hAx(7));
    bpLickRaster(SessionData, 8, 2, zeroField, '', hAx(8)); title('Uncued Big');    
    set(hAx, 'XLim', [lickHistSpecs(1), lickHistSpecs(2)]);     
    
%%     Photometry
    params.matpos = matpos_Ph;
    rows = 2;
    columns = 3;
    nAxes = 6;
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);

    

    

    
%     phPlotSessionAverage(analysis, type, outcome, fCh, varargin)    
    
    [~, hl] = phPlotSessionAverage(session.analysis, [1 2 3], [1 2 3], fCh, 'ax', hAx(1), 'linespec', {'g', 'k', 'r'});
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 2, fCh);
%     boundedline(avgX, avg, avgSEM, 'g', hAx(1), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 3, fCh);    
%     boundedline(avgX, avg, avgSEM, 'k', hAx(1), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 4, fCh);    
%     boundedline(avgX, avg, avgSEM, 'r', hAx(1), 'alpha');
    legend(hl, {'small', 'big', 'omit'}, 'Location', 'northwest');
    title('Short Tone');
    
    [~, hl] = phPlotSessionAverage(session.analysis, [4 5 6], [1 2 3], fCh, 'ax', hAx(2), 'linespec', {'g', 'k', 'r'});
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 2, fCh);
%     boundedline(avgX, avg, avgSEM, 'g', hAx(2), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 3, fCh);    
%     boundedline(avgX, avg, avgSEM, 'k', hAx(2), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 4, fCh);    
%     boundedline(avgX, avg, avgSEM, 'r', hAx(2), 'alpha');
    legend(hl, {'small', 'big', 'omit'}, 'Location', 'northwest');
    title('Long Tone');
    
    [~, hl] = phPlotSessionAverage(session.analysis, [1 4 7], [1 1 1], fCh, 'ax', hAx(3), 'linespec', {'c', 'm', 'b'});    
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 2, fCh);
%     boundedline(avgX, avg, avgSEM, 'c', hAx(3), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 2, fCh);    
%     boundedline(avgX, avg, avgSEM, 'm', hAx(3), 'alpha');
    legend(hl, {'short', 'long', 'uncued'}, 'Location', 'northwest');
    xlabel('time (s)'); title('small reward');    
    
    [~, hl] = phPlotSessionAverage(session.analysis, [2 5 8], [2 2 2], fCh, 'ax', hAx(4), 'linespec', {'c', 'm', 'b'});        
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 3, fCh);
%     boundedline(avgX, avg, avgSEM, 'c', hAx(4), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 3, fCh);    
%     boundedline(avgX, avg, avgSEM, 'm', hAx(4), 'alpha');
    legend(hl, {'short', 'long', 'uncued'}, 'Location', 'northwest');
    xlabel('time (s)'); title('big reward');    

    [~, hl] = phPlotSessionAverage(session.analysis, [3 6], [3 3], fCh, 'ax', hAx(5), 'linespec', {'c', 'm'});        
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 4, fCh);
%     boundedline(avgX, avg, avgSEM, 'c', hAx(5), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 4, fCh);    
%     boundedline(avgX, avg, avgSEM, 'm', hAx(5), 'alpha');
    legend(hl, {'small', 'big'}, 'Location', 'northwest');
    xlabel('time (s)'); title('omit');        
    
    
    
    
    
%% Saving
    saveas(fig, [fig_title '.fig']); %save as matlab fig
%     saveas(fig, [fig_title '.pdf']); %save as pdf
    disp('*** SOv2_Training_FS_SummaryFigure complete ***');
