
% using bpLoadSessions loaded the following data:
% ChAT_32, August 15 -> 26
% ChAT_31, August 15 -> 24
% ChAT_30, August 17 -> 24
% ChAT_29, August 17 -> 24
% ChAT_28, August 17 -> 23
% ChAT_26, August 17 -> 20, 22 -> 24
% and used makeTE_CuedOutcome_Odor_Complete to make a TE structure where each 
% element in the structure was a TE corresponding to a single mouse
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_CuedOutcome_Odor_Complete(sessions);
%%
TE.Photometry = processTrialAnalysis_Photometry2(sessions);
%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete_161007\ChAT_42';
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



%% generate trial lookups for different combinations of conditions

    highValueTrials = filterTE(TE, 'trialType', 1:3, 'reject', 0);
    lowValueTrials = filterTE(TE, 'trialType', 4:6, 'reject', 0);
    rewardTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    punishTrials = filterTE(TE, 'trialOutcome', 2, 'reject', 0);
    omitTrials = filterTE(TE, 'trialOutcome', 3, 'reject', 0);
    trialTypes = 1:9;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end
    %% plot photometry averages
    ensureFigure('Photometry_Averages', 1); 
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
    
    saveas(gcf, fullfile(savepath, 'phAverages_CuedOutcome_Odor_Complete_ChAT_26.fig'));   
    

    
    %% plot photometry rasters
    CLimFactor = 2;
    ensureFigure('phRastersFromTE_reward', 1); 
    subplot(1,4,1); phRasterFromTE(TE, trialsByType{1}, 1, 'CLimFactor', CLimFactor);
    title('hival, reward'); 
    subplot(1,4,2); phRasterFromTE(TE, trialsByType{4}, 1, 'CLimFactor', CLimFactor);
    title('loval, reward'); xlabel('time from reinforcement (s)'); 
    subplot(1,4,3); phRasterFromTE(TE, trialsByType{7}, 1, 'CLimFactor', CLimFactor);
    title('uncued, reward');
    subplot(1,4,4); phRasterFromTE(TE, trialsByType{3}, 1, 'CLimFactor', CLimFactor);
    title('hival, neutral'); xlabel('time from tone (s)'); 
    
    ensureFigure('phRastersFromTE_punish', 1); 
    subplot(1,4,1); phRasterFromTE(TE, trialsByType{2}, 1, 'CLimFactor', CLimFactor);
    title('hival, punish'); 
    subplot(1,4,2); phRasterFromTE(TE, trialsByType{5}, 1, 'CLimFactor', CLimFactor);
    title('loval, punish'); xlabel('time from reinforcement (s)'); 
    subplot(1,4,3); phRasterFromTE(TE, trialsByType{8}, 1, 'CLimFactor', CLimFactor);
    title('uncued, punish');
    subplot(1,4,4); phRasterFromTE(TE, trialsByType{6}, 1, 'CLimFactor', CLimFactor);
    title('loval, neutral'); xlabel('time from tone (s)'); 
   

    %% off conditions
    %% plot photometry rasters
    ensureFigure('phRasters_off', 1); 
    subplot(1,2,1); imshow(TE.Photometry.data(1).dFF(trialsByType{2}, :), [-.01 .01]); 
    colormap default; title('hival, punish'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2); imshow(TE.Photometry.data(1).dFF(trialsByType{4}, :), [-.01 .01]); 
    colormap default; title('loval, reward'); xlabel('time from reinforcement (s)');     

    %% plot photometry rasters 2 lab meeting
    ensureFigure('phRasters_hival', 1); 
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
    saveas(gcf, fullfile(savepath, 'phRasters_hiVal.fig'));  
    
    %% plot photometry rasters lab meeting low value
    ensureFigure('phRasters_lowVal', 1); 
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
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]);        
    set(gca, 'FontSize', 14)
    saveas(gcf, fullfile(savepath, 'phRasters_lowVal.fig'));  
    %% plot lick rasters
    ensureFigure('lickRasters', 1);
    subplot(1,2,1);
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]);
    
    subplot(1,2,2);
    eventRasterFromTE(TE, trialsByType{5}, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('lowval, punish'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]);
    
    
