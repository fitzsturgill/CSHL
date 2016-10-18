

%% 
sessions = bpLoadSessions;
%%
TE = makeTE_CuedOutcome_Odor_Complete(sessions);
%%
TE.Photometry = processTrialAnalysis_Photometry2(sessions);

%% extract peak trial dFF responses to cues and reinforcement
TE.phPeak_cs = bpCalcPeak_dFF(TE.Photometry, 1, [0 3], TE.Cue, 'method', 'peak');

%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
subjectName = TE.filename{1}(1:7);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);
%%
truncateSessionsFromTE(TE, 'init');
%%
save(fullfile(savepath, 'TE.mat'), 'TE');
disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
%% make tiled array of antic. licks for low and high value odors vs trial number
smoothFactor = 11;
ensureFigure('AnticLickRate_crossSessions', 1);

highTrials = filterTE(TE, 'trialType', 1:3, 'reject', 0);
lowTrials = filterTE(TE, 'trialType', 4:6, 'reject', 0);
rewardTrials = filterTE(TE, 'trialOutcome', 1);    
axes; plot(find(highTrials), smooth(TE.csLicks.rate(highTrials), smoothFactor), 'b.'); hold on; 
plot(find(lowTrials), smooth(TE.csLicks.rate(lowTrials), smoothFactor), 'r.');
plot(find(rewardTrials), smooth(TE.usLicks.rate(rewardTrials), smoothFactor), 'k.')    
plot(1:length(TE.trialNumber), [0; diff(TE.sessionIndex)] * 10, 'g');
ylabel('Lick/s, Us'); xlabel('trial #'); textBox(TE.filename{1}(1:7));
set(gca, 'YLim', [0 10]);

saveas(gcf, fullfile(savepath, 'AnticLickRate_crossSessions.fig'));
saveas(gcf, fullfile(savepath, 'AnticLickRate_crossSessions.jpg'));



%% generate trial lookups for different combinations of conditions

    highValueTrials = filterTE(TE, 'trialType', 1:3, 'reject', 0);
    lowValueTrials = filterTE(TE, 'trialType', 4:6, 'reject', 0);
    uncuedTrials = filterTE(TE, 'trialType', 7:9, 'reject', 0);    
    rewardTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    punishTrials = filterTE(TE, 'trialOutcome', 2, 'reject', 0);    
    omitTrials = filterTE(TE, 'trialOutcome', 3, 'reject', 0);
    trialTypes = 1:9;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end
    %% plot photometry averages
    h=ensureFigure('Photometry_Averages', 1); 
    mcLandscapeFigSetup(h);
    pm = [3 2];
    
    % - 6 0 4
    subplot(pm(1), pm(2), 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 3 7]), 1); %high value, reward
    legend(hl, {'hival, rew', 'hival, omit', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('high value'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); textBox(TE.filename{1}(1:7));

    subplot(pm(1), pm(2), 2); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 6 8]), 1); % low value, punish
    legend(hl, {'loval, pun', 'loval, omit', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('low value'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 3); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1); % reward, varying degrees of expectation
    legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('reward all'); ylabel('dF/F'); xlabel('time from reinforcement (s)');     

    subplot(pm(1), pm(2), 4); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 2 8]), 1); % punishment, varying degrees of expectation
    legend(hl, {'loval, pun', 'hival, pun', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('punish all'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 5); [ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
        'window', [-6 0], 'linespec', {'m', 'g'}); hold on;
    
    subplot(pm(1), pm(2), 5); [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, omitTrials}, 1,...
        'window', [0 4], 'linespec', {'b', 'r', 'k'});
    hl = [hla hl];
    legend(hl, {'loval', 'hival', 'rew', 'pun', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Balazs'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 
    
    
    saveas(gcf, fullfile(savepath, 'phAverages.fig'));
    saveas(gcf, fullfile(savepath, 'phAverages.jpg'));
    
    %% plot lick averages
    h = ensureFigure('Lick_Averages', 1);
%     mcPortraitFigSetup(h);
    mcLandscapeFigSetup(h);
    pm = [2 2];
    
    % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-4 0], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
    axh = [];
    subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 'Port1In', varargin{:});
    legend(hl, {'hival', 'loval', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Cue Licks'); ylabel('licks (s)'); xlabel('time from reinforcement (s)'); textBox(TE.filename{1}(1:7));  

    % window changed to US
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-1 4], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
    
    % reward
    axh(end + 1) = subplot(pm(1), pm(2), 2); [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 4 7]), 'Port1In', varargin{:});
    legend(hl, {'hival', 'loval', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Reward'); ylabel('licks (s)'); xlabel('time r (s)');
    
    % punish
    axh(end + 1) = subplot(pm(1), pm(2), 3); [ha, hl] = plotEventAverageFromTE(TE, trialsByType([2 5 8]), 'Port1In', varargin{:});
    legend(hl, {'hival', 'loval', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Punish'); ylabel('licks (s)'); xlabel('time r (s)');
    
    % neutral
    axh(end + 1) = subplot(pm(1), pm(2), 4); [ha, hl] = plotEventAverageFromTE(TE, trialsByType([3 6 9]), 'Port1In', varargin{:});
    legend(hl, {'hival', 'loval', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Neutral'); ylabel('licks (s)'); xlabel('time r (s)');
    
    sameYScale(axh) % match y scaling
    saveas(gcf, fullfile(savepath, 'lickAverages.fig'));
    saveas(gcf, fullfile(savepath, 'lickAverages.jpg'));    
   
  

    %% plot photometry rasters
    CLimFactor = 2;
    h=ensureFigure('phRastersFromTE_reward', 1);
    mcPortraitFigSetup(h);
    subplot(1,4,1); phRasterFromTE(TE, trialsByType{1}, 1, 'CLimFactor', CLimFactor);
    title([TE.filename{1}(1:7) ': hival, reward'], 'Interpreter', 'none'); 
    subplot(1,4,2); phRasterFromTE(TE, trialsByType{4}, 1, 'CLimFactor', CLimFactor);
    title('loval, reward'); xlabel('time from reinforcement (s)'); 
    subplot(1,4,3); phRasterFromTE(TE, trialsByType{7}, 1, 'CLimFactor', CLimFactor);
    title('uncued, reward');
    subplot(1,4,4); phRasterFromTE(TE, trialsByType{3}, 1, 'CLimFactor', CLimFactor);
    title('hival, neutral'); xlabel('time from tone (s)'); 

    saveas(gcf, fullfile(savepath, 'phRasters_reward.fig'));
    saveas(gcf, fullfile(savepath, 'phRasters_reward.jpg'));    
    
    h = ensureFigure('phRastersFromTE_punish', 1); 
    mcPortraitFigSetup(h);
    subplot(1,4,1); phRasterFromTE(TE, trialsByType{2}, 1, 'CLimFactor', CLimFactor);
    title([TE.filename{1}(1:7) ': hival, punish'], 'Interpreter', 'none'); 
    subplot(1,4,2); phRasterFromTE(TE, trialsByType{5}, 1, 'CLimFactor', CLimFactor);
    title('loval, punish'); xlabel('time from reinforcement (s)'); 
    subplot(1,4,3); phRasterFromTE(TE, trialsByType{8}, 1, 'CLimFactor', CLimFactor);
    title('uncued, punish');
    subplot(1,4,4); phRasterFromTE(TE, trialsByType{6}, 1, 'CLimFactor', CLimFactor);
    title('loval, neutral'); xlabel('time from tone (s)'); 
  
    saveas(gcf, fullfile(savepath, 'phRasters_punish.fig'));
    saveas(gcf, fullfile(savepath, 'phRasters_punish.jpg'));    

    %% plot photometry rasters 2 lab meeting
    h=ensureFigure('phRasters_hival', 1);
    mcPortraitFigSetup(h);
    prt = trialsByType{1};
    prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,2,1); 
    image('Xdata', [-6 4], 'YData', [1 find(length(prt))],...
        'CData', prcd, 'CDataMapping', 'Scaled', 'Parent', gca);
    set(gca, 'CLim', [-0.01 .01], 'YDir', 'Reverse');
    set(gca, 'XLim', [-6 4]);
    set(gca, 'FontSize', 14)
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2);
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]); 
    set(gca, 'FontSize', 14)
    saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_reward.fig'));
    saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_reward.jpg'));    
    
    %% plot photometry rasters lab meeting low value
    h=ensureFigure('phRasters_lowVal', 1); 
    mcPortraitFigSetup(h);    
    prt = trialsByType{5};
    prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,2,1); 
    image('Xdata', [-6 4], 'YData', [1 find(length(prt))],...
        'CData', prcd, 'CDataMapping', 'Scaled', 'Parent', gca);
    set(gca, 'CLim', [-0.01 .01], 'YDir', 'Reverse');
    set(gca, 'XLim', [-6 4]);
    set(gca, 'FontSize', 14)
    title('loVal, punish'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2);
    eventRasterFromTE(TE, prt, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('loVal, punish'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]);        
    set(gca, 'FontSize', 14)
    saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_punish.fig'));
    saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_punish.jpg'));    


    
