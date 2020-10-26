% SO_RewardPunish_Odor_TE_analysisScript

saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_SO_RewardPunish_Odor_v2(sessions);
%%
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'simple', 'blMode', 'byTrial', 'zeroField', 'PostUsRecording');

%%
basepath = uigetdir;
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

%%
TE.csLicks = countEventFromTE(TE, 'Port1In', [-2 0], TE.PostUsRecording);
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.PostUsRecording);
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'ReinforcementOutcome', 'Reward'));
%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% generate trial lookups for different combinations of conditions
    validTrials = filterTE(TE, 'reject', 0);
    rewardOdorTrials = filterTE(TE, 'trialType', 1:2, 'reject', 0);
    punishOdorTrials = filterTE(TE, 'trialType', 3:4, 'reject', 0);
    uncuedTrials = filterTE(TE, 'trialType', 5:7, 'reject', 0);    
    rewardTrials = filterTE(TE, 'trialOutcome', [1 5], 'reject', 0);
    punishTrials = filterTE(TE, 'trialOutcome', [3 6], 'reject', 0);    
    omitTrials = filterTE(TE, 'trialOutcome', [2 4 7], 'reject', 0);
    trialTypes = 1:7;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end

    
%% averages
    channel = 1;
    fdField = 'ZS';
    window = [-7 4];
    saveName = sprintf('Averages_%s', subjectName);  
    ensureFigure(saveName, 1);
    
    subplot(2, 2, 1);
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([7 1 2 5]), channel,...
        'PhotometryField', 'Photometry', 'FluorDataField', fdField, 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'k', 'b','g','c'}); 
    legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    title('Reward'); ylabel(sprintf('BF %s', fdField)); textBox(subjectName, gca, [0.25 0.95], 8);%set(gca, 'YLim', ylim);
    
    subplot(2, 2, 2);
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([7 3 4 6]), channel,...
        'PhotometryField', 'Photometry', 'FluorDataField', fdField, 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'k', 'r','g','m'}); 
    legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    title('Punish'); ylabel(sprintf('BF %s', fdField)); textBox(subjectName, gca, [0.25 0.95], 8);%set(gca, 'YLim', ylim);
    
    subplot(2,2,3);
    varargin = {'window', [-4 0], 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'k', 'b','g', 'c'}};
    axh = [];
    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([7 1 2 5]), 'Port1In', varargin{:});
    legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    ylabel('licks (1/s)'); xlabel('time from reinforcement (s)'); 
    
    subplot(2,2,4);
    varargin = {'window', [-4 0], 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'k', 'r','g','m'}};
    axh = [];
    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([7 3 4 6]), 'Port1In', varargin{:});
    legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    ylabel('licks (1/s)'); xlabel('time from reinforcement (s)');   
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end
%%
    %% raster plots
    saveName = sprintf('rasters_%s', subjectName); 
    ensureFigure(saveName, 1);
    mcLandscapeFigSetup(gcf); 
    photometryField = 'Photometry';
    climfactor = 4;
    trialTypes = [1 2 5 3 4 6];
    titles  = {'cued rew.', 'omit rew.', 'uncued rew.', 'cued pun.', 'omit pun.', 'uncued pun.'};
    for typeCounter = 1:length(trialTypes)
        trialType = trialTypes(typeCounter);
        sessionChanges = find(diff(TE.sessionIndex(trialsByType{trialType}, :))) + 1;
        subplot(2,6,typeCounter);
        phRasterFromTE(TE, trialsByType{trialType}, channel, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor,...
            'zeroTimes', TE.usZeros', 'window', window, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,       
        title(titles{typeCounter}); 
        
        
        subplot(2,6,typeCounter + 6);
        eventRasterFromTE(TE, trialsByType{trialType}, 'Port1In', 'trialNumbering', 'consecutive',...
            'zeroTimes', TE.usZeros, 'window', window);
        line(repmat([window(1); window(2)], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
        set(gca, 'XLim', window); 
    end
    subplot(2,6,8); xlabel('Time from Reinforcement (s)');
    subplot(2,6,1); t = textBox(subjectName, gca, [0.1 0.95]); set(t, 'Color', [1 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end
%%