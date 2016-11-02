
saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_CuedOutcome_Odor_Complete(sessions);
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial');

%% extract peak trial dFF responses to cues and reinforcement and lick counts
TE.phPeak_cs = bpCalcPeak_dFF(TE.Photometry, 1, [0 2], TE.Cue, 'method', 'mean');
TE.phPeak_us = bpCalcPeak_dFF(TE.Photometry, 1, [0 0.5], TE.Us, 'method', 'mean');

TE.csLicks = countEventFromTE(TE, 'Port1In', [-2 0], TE.Us);
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.Us);

%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
subjectName = TE.filename{1}(1:7);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);
%%
truncateSessionsFromTE(TE, 'init');
%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% cross sessions bleaching curve and dual exponential fits
ensureFigure('sessionBleach_Correction', 1);
plot(TE.Photometry.data(1).blF_raw, 'k'); hold on;
plot(TE.Photometry.data(1).blF, 'r');
if saveOn
    saveas(gcf, fullfile(savepath, 'sessionBleach_correction.fig'));
    saveas(gcf, fullfile(savepath, 'sessionBleach_correction.jpg'));
end
%% cross trial bleaching fits for each session plotted as axis array
ensureFigure('trialBleach_Correction', 1);
nSessions = length(TE.Photometry.bleachFit);
subA = ceil(sqrt(nSessions));
for counter = 1:nSessions
    subplot(subA, subA, counter);
    plot(TE.Photometry.bleachFit(counter).trialTemplate, 'k'); hold on;
    plot(TE.Photometry.bleachFit(counter).trialFit, 'r');
%     title(num2str(counter));    
end
if saveOn
    saveas(gcf, fullfile(savepath, 'trialBleach_Correction.fig'));
    saveas(gcf, fullfile(savepath, 'trialBleach_Correction.jpg'));
end

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
if saveOn
    saveas(gcf, fullfile(savepath, 'AnticLickRate_crossSessions.fig'));
    saveas(gcf, fullfile(savepath, 'AnticLickRate_crossSessions.jpg'));
end


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
    subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 3 7]), 1); %high value, reward
    legend(hl, {'hival, rew', 'hival, omit', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('high value'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); textBox(TE.filename{1}(1:7));

    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 6 8]), 1); % low value, punish
    legend(hl, {'loval, pun', 'loval, omit', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('low value'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1); % reward, varying degrees of expectation
    legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('reward all'); ylabel('dF/F'); xlabel('time from reinforcement (s)');     

    subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 2 8]), 1); % punishment, varying degrees of expectation
    legend(hl, {'loval, pun', 'hival, pun', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('punish all'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 5, 'FontSize', 12, 'LineWidth', 1); [ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
        'window', [-6 0], 'linespec', {'m', 'g'}); hold on;
    
    subplot(pm(1), pm(2), 5); [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, omitTrials}, 1,...
        'window', [0 4], 'linespec', {'b', 'r', 'k'});
    hl = [hla hl];
    legend(hl, {'loval', 'hival', 'rew', 'pun', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Balazs'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 
    
if saveOn    
    saveas(gcf, fullfile(savepath, 'phAverages.fig'));
    saveas(gcf, fullfile(savepath, 'phAverages.jpg'));
end
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
if saveOn    
    saveas(gcf, fullfile(savepath, 'lickAverages.fig'));
    saveas(gcf, fullfile(savepath, 'lickAverages.jpg'));    
end
  

    %% plot photometry rasters
    CLimFactor = 4;
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
if saveOn
    saveas(gcf, fullfile(savepath, 'phRasters_punish.fig'));
    saveas(gcf, fullfile(savepath, 'phRasters_punish.jpg'));    
end

%%
    %% plot photometry rasters reward- alternate for lab meeting
    CLimFactor = 2;
    h=ensureFigure('phRastersFromTE_reward', 1);
    mcPortraitFigSetup(h);
    

%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,4,1); 
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('hival, reward'); xlabel('time from reinforcement (s)'); ylabel('trial number');
    set(gca, 'XLim', [-6 4]); 
    set(gca, 'YLim', [0 length(find(trialsByType{1}))]);
    set(gca, 'FontSize', 14)

    
    subplot(1,4,2); phRasterFromTE(TE, trialsByType{1}, 1, 'CLimFactor', CLimFactor);
    title('hival, reward', 'Interpreter', 'none'); 
        set(gca, 'FontSize', 14)
    subplot(1,4,3); phRasterFromTE(TE, trialsByType{4}, 1, 'CLimFactor', CLimFactor);
    title('loval, reward'); 
        set(gca, 'FontSize', 14)
    subplot(1,4,4); phRasterFromTE(TE, trialsByType{7}, 1, 'CLimFactor', CLimFactor);
    title('uncued, reward');
        set(gca, 'FontSize', 14)


    saveas(gcf, fullfile(savepath, 'phRasters_reward_alternate.fig'));
    saveas(gcf, fullfile(savepath, 'phRasters_reward_alternate.jpg'));    

%% plot photometry rasters 2 lab meeting
    h=ensureFigure('phRasters_hival', 1);
    mcPortraitFigSetup(h);
    prt = trialsByType{1};
%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,2,1); 
    phRasterFromTE(TE, trialsByType{1}, 1, 'CLimFactor', CLimFactor);
%     image('Xdata', [-6 4], 'YData', [1 find(length(prt))],...
%         'CData', prcd, 'CDataMapping', 'Scaled', 'Parent', gca);
%     set(gca, 'CLim', [-0.01 .01], 'YDir', 'Reverse');
    set(gca, 'XLim', [-6 4]);
    set(gca, 'YLim', [0 length(find(trialsByType{1}))]);
    set(gca, 'FontSize', 14)
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2);
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]); 
    set(gca, 'YLim', [0 length(find(trialsByType{1}))]);
    set(gca, 'FontSize', 14)
if saveOn
    saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_reward.fig'));
    saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_reward.jpg'));    
end
    %% plot photometry rasters lab meeting low value
    h=ensureFigure('phRasters_lowVal', 1); 
    mcPortraitFigSetup(h);    
    prt = trialsByType{5};
%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,2,1); 
    phRasterFromTE(TE, trialsByType{5}, 1, 'CLimFactor', CLimFactor);    
%     image('Xdata', [-6 4], 'YData', [1 find(length(prt))],...
%         'CData', prcd, 'CDataMapping', 'Scaled', 'Parent', gca);
    set(gca, 'CLim', [-0.01 .01], 'YDir', 'Reverse');
    set(gca, 'YLim', [0 length(find(trialsByType{5}))]);   
    set(gca, 'XLim', [-6 4]);
    set(gca, 'FontSize', 14)
    title('loVal, punish'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2);
    eventRasterFromTE(TE, prt, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('loVal, punish'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]);     
    set(gca, 'YLim', [0 length(find(trialsByType{5}))]);        
    set(gca, 'FontSize', 14)
if saveOn
    saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_punish.fig'));
    saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_punish.jpg'));    
end
    %% summary statistics
    cComplete_summary = struct(...
        'phCue_CV', zeros(1,9),...
        'phCue_avg', zeros(1,9),...
        'phOutcome_CV', zeros(1,9),...
        'phOutcome_avg', zeros(1,9),...
        'cueLicks_low', 0,...
        'cueLicks_high', 0,...
        'rewardLicks', 0);
    for counter = 1:length(trialTypes)
        trials = trialsByType{counter};
        peaks = TE.phPeak_cs.data(trials);
        cComplete_summary.phCue_CV(counter) = nanmean(peaks) / std(peaks, 'omitnan');
        cComplete_summary.phCue_avg(counter) = nanmean(peaks);
        peaks = TE.phPeak_us.data(trials);
        cComplete_summary.phOutcome_CV(counter) = nanmean(peaks) / std(peaks, 'omitnan');
        cComplete_summary.phOutcome_avg(counter) = nanmean(peaks); 
    end
    cComplete_summary.cueLicks_low = nanmean(TE.csLicks.rate(lowValueTrials));
    cComplete_summary.cueLicks_high = nanmean(TE.csLicks.rate(highValueTrials));        
    cComplete_summary.rewardLicks = nanmean(TE.usLicks.rate(rewardTrials));    
if saveOn
    save(fullfile(savepath, ['summary_' subjectName '.mat']), 'cComplete_summary');
    disp(['*** saving: ' fullfile(savepath, ['summary_' subjectName '.mat']) ' ***']);
end

    
    

    
    %% dFF vs licks scatter plot
    ensureFigure('phVSLicks_Scatter', 1);
    scatter(TE.usLicks.count(trialsByType{1}) + rand(length(find(trialsByType{1})), 1) - 0.5, TE.phPeak_us.data(trialsByType{1}), 'b'); hold on;
    scatter(TE.usLicks.count(trialsByType{4}) + rand(length(find(trialsByType{4})), 1) - 0.5, TE.phPeak_us.data(trialsByType{4}), 'r');
        set(gca, 'FontSize', 12); xlabel('reward licks (jittered)'); ylabel('phReward (dFF-avg)');
    set(gca, 'XLim', [0 30]);
if saveOn
    saveas(gcf, fullfile(savepath, 'dFF_vs_licks.fig'));
    saveas(gcf, fullfile(savepath, 'dFF_vs_licks.jpg'));        
end
    %% us vs cs dFF for highValue reward condition
    ensureFigure('phUSvsCS_Scatter', 1);
    scatter(TE.phPeak_cs.data(trialsByType{1}), TE.phPeak_us.data(trialsByType{1}), 'b'); hold on;
if saveOn
    saveas(gcf, fullfile(savepath, 'us_vs_cs_dFF.fig'));
    saveas(gcf, fullfile(savepath, 'us_vs_cs_dFF.jpg'));    
end
    %% is the trough related to number of licks? 
    TE.phTrough_us = bpCalcPeak_dFF(TE.Photometry, 1, [0 2], TE.Us, 'method', 'min');
    ensureFigure('phDip_Scatter', 1);
    scatter(TE.usLicks.count(rewardTrials) + rand(length(find(rewardTrials)), 1) - 0.5, TE.phTrough_us.data(rewardTrials), 'b'); 
    
    %% 
    


    
