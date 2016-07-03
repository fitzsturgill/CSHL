function phSummaryFigure_noPunish_v1(session)
% make and save summary figure for a session

% assumes that you've loaded session using bpLoadSession which creates
% session folder and CDs into that folder
    zeroField = 'DeliverStimulus';
    startField = 'ITI';
    endField = 'PostTrialRecording';
    endField = 'ResetDrinkingTimer';
    lickHistSpecs = [-1, 4, 0.1]; % start, stop and bin size of lick histograms
    fig = ensureFigure('phSummaryFigure_noPunish_v1', 1); % ensure figure, erase if preexisting
    fig=mcPortraitFigSetup(fig);


    
%% Define the axes matrix positions on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    matpos_title = [0 0 1 .1];
    matpos_rasters = [0 .1 1 .4];
    matpos_hist = [0 .5 1 .5];
    params.cellmargin = [.05 .05 0.05 0.05];    
    
%% Figure Title
    [~, fig_title, ~] = fileparts(session.filename); % session name
    title_ax = textAxes(fig, fig_title, 12);
    params.matpos = matpos_title;
    setaxesOnaxesmatrix(title_ax, 1, 1, 1, params, fig);
    
%% Lick Rasters
%     trialTypes = [1 2; 1 3; 1 4; 2 2; 2 3; 2 4];
    params.matpos = matpos_rasters;
    rows = 2;
    columns = 2;
    nAxes = 4;
    % layout axes in grid
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);
    bpLickRaster(session.SessionData, 1, 2, zeroField, '', hAx(1)); ylabel('Tone A'); title('Reward');
    bpLickRaster(session.SessionData, 1, 3, zeroField, '', hAx(2)); title('Omit');
%     bpLickRaster(session.SessionData, 1, 4, zeroField, '', hAx(3)); title('Punish');
    bpLickRaster(session.SessionData, 2, 2, zeroField, '', hAx(3)); ylabel('Tone B');
    bpLickRaster(session.SessionData, 2, 3, zeroField, '', hAx(4)); 
%     bpLickRaster(session.SessionData, 2, 4, zeroField, '', hAx(6));
    
    set(hAx, 'XLim', [lickHistSpecs(1), lickHistSpecs(2)]);

%% Histograms    
    
    params.matpos = matpos_hist;
    rows = 1;
    columns = 2;
    nAxes = 2;
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);
    % reward vs punishment lick histogram
    bpLickHist(session.SessionData, [1 2], [2 2], lickHistSpecs, zeroField, startField, endField, {'g', 'r'}, '', hAx(1));
    xlabel('Time (s) from tone'); ylabel('Lick Rate');
%     % reward vs punishment photometry data
% %     phPlotSessionAvg(session.SessionData, [1 2], [2 4], '', hAx(2), {'g', 'r'});
% %     xlabel('Time (s) from tone'); ylabel('dF/F');


%     phPlotSessionAvg(session.SessionData, [1 1], [2 3], '', hAx(1), {'g', 'k', 'r'});
%     xlabel('Time (s) from tone'); ylabel('dF/F'); textBox('Tone A');
%     phPlotSessionAvg(session.SessionData, [2 2], [2 3], '', hAx(2), {'g', 'k', 'r'});
%     xlabel('Time (s) from tone'); ylabel('dF/F'); textBox('Tone B');
    
    
%% Saving
    saveas(fig, [fig_title '.fig']); %save as matlab fig
    saveas(fig, [fig_title '.pdf']); %save as pdf
    disp('*** phSummaryFigure_noPunish_v1 complete ***');
    
    
    