function session = SOv2_Training_FS_SummaryFigure(session)
% make and save summary figure for the SOv2_Training_FS protocol

% I haven't yet updated this to use processAnalysis_Event

% assumes that 1) you've loaded session using session = bpLoadSession 


    zeroField = 'DeliverStimulus';
    startField = 'ITI';
%     if noFluor
%         endField = {'PostUS', 'last'};
%     else
        endField = 'PostTrialRecording';
%     end
    lickHistSpecs = [-2, 4, 0.25]; % start, stop and bin size of lick histograms
    fig = ensureFigure('SOv2_Training_FS_SummaryFigure', 1); % ensure figure, erase if preexisting
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
    bpLickRaster(SessionData, 1, 2, zeroField, '', hAx(1)); ylabel('Tone A'); title('Reward');
    bpLickRaster(SessionData, 1, 3, zeroField, '', hAx(2)); title('Omit');
    bpLickRaster(SessionData, 1, 4, zeroField, '', hAx(3)); title('Punish');
    bpLickRaster(SessionData, 3, 2, zeroField, '', hAx(4)); title('Uncued Reward');    
    bpLickRaster(SessionData, 2, 2, zeroField, '', hAx(5)); ylabel('Tone B'); xlabel('time (s) from tone');
    bpLickRaster(SessionData, 2, 3, zeroField, '', hAx(6)); 
    bpLickRaster(SessionData, 2, 4, zeroField, '', hAx(7));
    bpLickRaster(SessionData, 4, 4, zeroField, '', hAx(8)); title('Uncued Punish');    
    set(hAx, 'XLim', [lickHistSpecs(1), lickHistSpecs(2)]);    
    
%     % histograms
%     % tone A + reward, punishment and ommission 
%     trialTypes = [1 1 1];
%     trialOutcomes = {[1 2] [0 4] [3 5]};    
%     bpLickHist(SessionData, trialTypes, trialOutcomes, lickHistSpecs, zeroField, startField, endField, {'g', 'r', 'k'}, '', hAx(4));
%     title('Tone A'); xlabel('Time (s) from tone');
%     
%     % tone B + reward, punishment and ommission 
%     trialTypes = [2 2 2];
%     trialOutcomes = {[1 2] [0 4] [3 5]};   
%     bpLickHist(SessionData, trialTypes, trialOutcomes, lickHistSpecs, zeroField, startField, endField, {'g', 'r', 'k'}, '', hAx(8));
%     title('Tone B');xlabel('Time (s) from tone'); ylabel('Lick Rate');        
%     set(hAx, 'XLim', [lickHistSpecs(1), lickHistSpecs(2)]);

%% Photometry

    params.matpos = matpos_Ph;
    rows = 2;
    columns = 3;
    nAxes = 6;
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);

    
    fCh = 1;
    
    conditions = {'toneA_reward', 'toneA_omit', 'toneA_punish', 'toneB_reward', 'toneB_omit', 'toneB_punish'};
    trialTypes = [ 1 1 1 2 2 2];
    outcomes = [2 3 4 2 3 4];
    linespecs = {'g', 'k', 'r', 'g', 'k', 'r'};
    
%     phPlotSessionAverage(analysis, type, outcome, fCh, varargin)    
    
    [~, hl] = phPlotSessionAverage(session.analysis, [1 1 1], [2 3 4], fCh, 'ax', hAx(1), 'linespec', {'g', 'k', 'r'});
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 2, fCh);
%     boundedline(avgX, avg, avgSEM, 'g', hAx(1), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 3, fCh);    
%     boundedline(avgX, avg, avgSEM, 'k', hAx(1), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 4, fCh);    
%     boundedline(avgX, avg, avgSEM, 'r', hAx(1), 'alpha');
    legend(hl, {'reward', 'omit', 'punish'}, 'Location', 'northwest');
    title('tone A');
    
    [~, hl] = phPlotSessionAverage(session.analysis, [2 2 2], [2 3 4], fCh, 'ax', hAx(2), 'linespec', {'g', 'k', 'r'});
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 2, fCh);
%     boundedline(avgX, avg, avgSEM, 'g', hAx(2), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 3, fCh);    
%     boundedline(avgX, avg, avgSEM, 'k', hAx(2), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 4, fCh);    
%     boundedline(avgX, avg, avgSEM, 'r', hAx(2), 'alpha');
    legend(hl, {'reward', 'omit', 'punish'}, 'Location', 'northwest');
    title('tone B');
    
    [~, hl] = phPlotSessionAverage(session.analysis, [1 2 3], [2 2 2], fCh, 'ax', hAx(3), 'linespec', {'c', 'm', 'b'});    
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 2, fCh);
%     boundedline(avgX, avg, avgSEM, 'c', hAx(3), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 2, fCh);    
%     boundedline(avgX, avg, avgSEM, 'm', hAx(3), 'alpha');
    legend(hl, {'toneA', 'toneB', 'uncued'}, 'Location', 'northwest');
    xlabel('time (s)'); title('reward');    
    
    [~, hl] = phPlotSessionAverage(session.analysis, [1 2], [3 3], fCh, 'ax', hAx(4), 'linespec', {'c', 'm'});        
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 3, fCh);
%     boundedline(avgX, avg, avgSEM, 'c', hAx(4), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 3, fCh);    
%     boundedline(avgX, avg, avgSEM, 'm', hAx(4), 'alpha');
    legend(hl, {'toneA', 'toneB'}, 'Location', 'northwest');
    xlabel('time (s)'); title('omit');    

    [~, hl] = phPlotSessionAverage(session.analysis, [1 2 4], [4 4 4], fCh, 'ax', hAx(5), 'linespec', {'c', 'm', 'b'});        
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 1, 4, fCh);
%     boundedline(avgX, avg, avgSEM, 'c', hAx(5), 'alpha');
%     [avg, avgX, avgSEM] = phSessionAverage(SessionData, 2, 4, fCh);    
%     boundedline(avgX, avg, avgSEM, 'm', hAx(5), 'alpha');
    legend(hl, {'toneA', 'toneB', 'uncued'}, 'Location', 'northwest');
    xlabel('time (s)'); title('punish');        
    
    
    
    
    
%% Saving
    saveas(fig, [fig_title '.fig']); %save as matlab fig
%     saveas(fig, [fig_title '.pdf']); %save as pdf
    disp('*** SOv2_Training_FS_SummaryFigure complete ***');
    
    
    