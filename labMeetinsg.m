
% make and save summary figure for a session

% assumes that you've loaded session using bpLoadSession which creates
% session folder and CDs into that folder


    zeroField = 'DeliverStimulus';
    startField = 'ITI';
    endField = 'PostUS';
%     endField = {'PostUS', 'last'};

    lickHistSpecs = [-1, 4, 0.1]; % start, stop and bin size of lick histograms
    fig = ensureFigure('phSummaryFigure_v1', 1); % ensure figure, erase if preexisting
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
    columns = 3;
    nAxes = 6;
    % layout axes in grid
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);
    bpLickRaster(session.SessionData, 1, [1 2], zeroField, '', hAx(1)); ylabel('Tone A'); title('Reward');
    bpLickRaster(session.SessionData, 1, [0 4], zeroField, '', hAx(2)); title('Punish');
    bpLickRaster(session.SessionData, 1, [3 5], zeroField, '', hAx(3)); title('Omit');
    bpLickRaster(session.SessionData, 2, [1 2], zeroField, '', hAx(4)); ylabel('Tone B');
    bpLickRaster(session.SessionData, 2, [0 4], zeroField, '', hAx(5));
    bpLickRaster(session.SessionData, 2, [3 5], zeroField, '', hAx(6))
    
    set(hAx, 'XLim', [lickHistSpecs(1), lickHistSpecs(2)]);

%% Histograms    
    
    params.matpos = matpos_hist;
    rows = 2;
    columns = 2;
    nAxes = 4;
    hAx = axesmatrix(rows, columns, 1:nAxes, params, fig);
    trialTypes = [1 2];
    trialOutcomes = {[1 2] [0 4]};
%     % tone A + reward vs tone B + punishment lick histogram
%     bpLickHist(session.SessionData, trialTypes, trialOutcomes, lickHistSpecs, zeroField, startField, endField, {'g', 'r'}, '', hAx(1));
%     title('Tone A vs B'); ylabel('Lick Rate');
    
    % tone A + reward, punishment and ommission 
    trialTypes = [1 1 1];
    trialOutcomes = {[1 2] [0 4] [3 5]};    
    bpLickHist(session.SessionData, trialTypes, trialOutcomes, lickHistSpecs, zeroField, startField, endField, {'g', 'r', 'k'}, '', hAx(1));
    title('Tone A'); xlabel('Time (s) from tone');
    
    % tone B + reward, punishment and ommission 
    trialTypes = [2 2 2];
    trialOutcomes = {[1 2] [0 4] [3 5]};   
    bpLickHist(session.SessionData, trialTypes, trialOutcomes, lickHistSpecs, zeroField, startField, endField, {'g', 'r', 'k'}, '', hAx(2));
    title('Tone B');xlabel('Time (s) from tone'); ylabel('Lick Rate');    
    
    
% plot fluo data
conditions = {'toneA_reward', 'toneA_omit', 'toneA_punish', 'toneB_reward', 'toneB_omit', 'toneB_punish'};
trialTypes = [ 1 1 1 2 2 2];
outcomes = [2 3 4 2 3 4];
linespecs = {'g', 'k', 'r', 'g', 'k', 'r'};

thisAx = hAx(3);

for counter = 1:length(conditions)


    if counter == 4
      thisAx = hAx(4);
    end
    boundedline(demodData(counter).avgX, demodData(counter).avg, demodData(counter).avgSEM, linespecs{counter}, gca, 'alpha');     
    if counter == 1
        ylabel('GCAMP6s'); xlabel('time(s)'); title('tone A');        
    elseif counter == 4
        ylabel('GCAMP6s'); xlabel('time(s)'); title('tone B');        
    end
end
    
    
%% Saving
    saveas(fig, [fig_title 'labMeeting' '.fig']); %save as matlab fig
    saveas(fig, [fig_title 'labMeeting' '.pdf']); %save as pdf
    disp('*** phSummaryFigure_v1 complete ***');
    
    
    