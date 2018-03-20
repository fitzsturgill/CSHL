
saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
tau = [1 1];
%%
TE = makeTE_CuedOutcome_Odor_Complete(sessions);
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial', 'tau', tau);

%% if you want to try baseline expFit option... 
TE = makeTE_CuedOutcome_Odor_Complete(sessions);
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'simple', 'blMode', 'expFit');
%% For Dopamine GCaMP recordings, don't use expfit dFFMode option
TE = makeTE_CuedOutcome_Odor_Complete(sessions);
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'simple', 'blMode', 'byTrial');

%% extract peak trial dFF responses to cues and reinforcement and lick counts
TE.phPeak_cs = bpCalcPeak_dFF(TE.Photometry, 1, [0 2], TE.Cue, 'method', 'mean', 'phField', 'ZS');
TE.phPeak_us = bpCalcPeak_dFF(TE.Photometry, 1, [0 0.5], TE.Us, 'method', 'mean', 'phField', 'ZS');
TE.phPeak_preUs = bpCalcPeak_dFF(TE.Photometry, 1, [-0.5 0], TE.Us, 'method', 'mean', 'phField', 'ZS'); % to serve as local baseline to show positive reward and punishment responses
TE.phPeak_us_med = bpCalcPeak_dFF(TE.Photometry, 1, [0 0.5], TE.Us, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS'); % try median too
TE.phPeak_preUs_med = bpCalcPeak_dFF(TE.Photometry, 1, [-0.5 0], TE.Us, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS'); % median,  % to serve as local baseline to show positive reward and punishment responses
TE.phPeak_cs_phasic = bpCalcPeak_dFF(TE.Photometry, 1, [0 1], TE.Cue, 'method', 'mean', 'phField', 'ZS');
TE.phPeak_cs_sustained = bpCalcPeak_dFF(TE.Photometry, 1, [1 3], TE.Cue, 'method', 'mean', 'phField', 'ZS');

TE.csLicks = countEventFromTE(TE, 'Port1In', [-2 0], TE.Us);
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.Us);


%%
basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
if length(unique(TE.filename)) > 1
    sep = strfind(TE.filename{1}, '_');
    subjectName = TE.filename{1}(1:sep(2)-1);
else
    sep = strfind(TE.filename{1}, '.');
    subjectName = TE.filename{1}(1:sep(1)-1);
end
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);
%% truncate sessions according to when mouse becomes sated
rewardTrialsTrunc = filterTE(TE, 'trialOutcome', 1);
truncateSessionsFromTE(TE, 'init', 'usLicks', rewardTrialsTrunc);

%%
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', 11, 'Fs', 20, 'startField', 'Start');
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
    validTrials = filterTE(TE, 'reject', 0);
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
    % plot photometry averages
    h=ensureFigure('Photometry_Averages', 1); 
    mcLandscapeFigSetup(h);

    pm = [3 2];
    
    % - 6 0 4
    subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 3 7]), 1, 'FluorDataField', 'ZS'); %high value, reward
    legend(hl, {'hival, rew', 'hival, omit', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('high value'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); textBox(TE.filename{1}(1:7));

    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 6 8]), 1, 'FluorDataField', 'ZS'); % low value, punish
    legend(hl, {'loval, pun', 'loval, omit', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('low value'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1, 'FluorDataField', 'ZS'); % reward, varying degrees of expectation
    legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('reward all'); ylabel('Z Score'); xlabel('time from reinforcement (s)');     

    subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 2 8]), 1, 'FluorDataField', 'ZS'); % punishment, varying degrees of expectation
    legend(hl, {'loval, pun', 'hival, pun', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('punish all'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 5, 'FontSize', 12, 'LineWidth', 1); [ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
        'window', [-6 0], 'linespec', {'m', 'g'}, 'FluorDataField', 'ZS'); hold on;
    
    subplot(pm(1), pm(2), 5); [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, omitTrials}, 1,...
        'window', [0 4], 'linespec', {'b', 'r', 'k'}, 'FluorDataField', 'ZS');
    hl = [hla hl];
    legend(hl, {'loval', 'hival', 'rew', 'pun', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Balazs'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 
    
    subplot(pm(1), pm(2), 6, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([3 6 9]), 1, 'FluorDataField', 'ZS'); % reward, varying degrees of expectation
    legend(hl, {'hival, neutral', 'loval, neutral', 'neutral'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('neutral all'); ylabel('Z Score'); xlabel('time from reinforcement (s)');    
    
if saveOn    
    saveas(gcf, fullfile(savepath, 'phAverages.fig'));
    saveas(gcf, fullfile(savepath, 'phAverages.jpg'));
end

%% graph to pick phasic and sustained analysis windows
ensureFigure('windowPick', 1);
axes;
phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
        'window', [-4 0], 'linespec', {'m', 'g'}, 'FluorDataField', 'ZS'); 
    

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
if saveOn
    saveas(gcf, fullfile(savepath, 'phRasters_reward.fig'));
    saveas(gcf, fullfile(savepath, 'phRasters_reward.jpg'));    
end
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
%     TE.phPeak_cs_phasic
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

%% summary stastics re-do- for CCN talk
%     TE.phPeak_cs_phasic
    values = struct(...
        'avg', [],...
        'n', [],...
        'std', [],...
        'sem', []...
        );
            
    cComplete_summary2 = struct(...
        'phCue_phasic_low', values,...
        'phCue_phasic_high', values,...
        'phCue_sustained_low', values,...
        'phCue_sustained_high', values,...        
        'cueLicks_low', values,...
        'cueLicks_high', values,...
        'phReward_mean', values,...
        'phReward_mean_bl', values,...
        'phReward_med', values,...
        'phReward_med_bl', values,...        
        'phPunish_mean', values,...
        'phPunish_mean_bl', values,...
        'phPunish_med', values,...
        'phPunish_med_bl', values...                
    );
    
    summary2_data = {...
        'phCue_phasic_low', TE.phPeak_cs_phasic.data(lowValueTrials);...
        'phCue_phasic_high', TE.phPeak_cs_phasic.data(highValueTrials);...
        'phCue_sustained_low', TE.phPeak_cs_sustained.data(lowValueTrials);...
        'phCue_sustained_high', TE.phPeak_cs_sustained.data(highValueTrials);...        
        'cueLicks_low', TE.csLicks.rate(lowValueTrials);...
        'cueLicks_high', TE.csLicks.rate(highValueTrials);...
        'phReward_mean', TE.phPeak_us.data(rewardTrials & uncuedTrials);...
        'phReward_mean_bl', TE.phPeak_us.data(rewardTrials & uncuedTrials) - TE.phPeak_preUs.data(rewardTrials & uncuedTrials);...
        'phReward_med', TE.phPeak_us_med.data(rewardTrials & uncuedTrials);...
        'phReward_med_bl', TE.phPeak_us_med.data(rewardTrials & uncuedTrials) - TE.phPeak_preUs_med.data(rewardTrials & uncuedTrials);...
        'phPunish_mean', TE.phPeak_us.data(punishTrials & uncuedTrials);...
        'phPunish_mean_bl', TE.phPeak_us.data(punishTrials & uncuedTrials) - TE.phPeak_preUs.data(punishTrials & uncuedTrials);...
        'phPunish_med', TE.phPeak_us_med.data(punishTrials & uncuedTrials);...
        'phPunish_med_bl', TE.phPeak_us_med.data(punishTrials & uncuedTrials) - TE.phPeak_preUs_med.data(punishTrials & uncuedTrials);...
    };

% 'avg', 'n', 'std', 'sem'
    for counter = 1:size(summary2_data, 1)
        peaks = summary2_data{counter, 2};
        field = summary2_data{counter, 1};
        cComplete_summary2.(field).avg = nanmean(peaks);        
        cComplete_summary2.(field).n = sum(~isnan(peaks));
        cComplete_summary2.(field).std = std(peaks, 'omitnan');         
        cComplete_summary2.(field).sem = std(peaks, 'omitnan') / sqrt(sum(~isnan(peaks)));
    end

if saveOn
    save(fullfile(savepath, ['summary2_' subjectName '.mat']), 'cComplete_summary2');
    disp(['*** saving: ' fullfile(savepath, ['summary2_' subjectName '.mat']) ' ***']);
end
    

%% summary stastics, outcome responses, all conditions, baselined
%     TE.phPeak_cs_phasic
    data = struct(...
        'avg', [],...
        'n', [],...
        'std', [],...
        'sem', [],...
        'avg_med', [],...
        'n_med', [],...
        'std_med', [],...
        'sem_med', []...        
        );
            
    outcome_bl_all = TE.phPeak_us.data - TE.phPeak_preUs.data; 
    outcome_med_bl_all = TE.phPeak_us_med.data - TE.phPeak_preUs_med.data;

    ccomplete_outcomeSummary = repmat(data, 1, 9);

    for type = 1:9
            trials = trialsByType{type};
            peaks = outcome_bl_all(trials);
            peaks_med = outcome_med_bl_all(trials);
            ccomplete_outcomeSummary(type).avg = nanmean(peaks);        
            ccomplete_outcomeSummary(type).n = sum(~isnan(peaks));
            ccomplete_outcomeSummary(type).std = std(peaks, 'omitnan');         
            ccomplete_outcomeSummary(type).sem = std(peaks, 'omitnan') / sqrt(sum(~isnan(peaks)));           
            ccomplete_outcomeSummary(type).avg_med = nanmean(peaks_med);        
            ccomplete_outcomeSummary(type).n_med = sum(~isnan(peaks_med));
            ccomplete_outcomeSummary(type).std_med = std(peaks_med, 'omitnan');         
            ccomplete_outcomeSummary(type).sem_med = std(peaks_med, 'omitnan') / sqrt(sum(~isnan(peaks_med)));                       
    end

if saveOn
    save(fullfile(savepath, ['ccomplete_outcomeSummary_' subjectName '.mat']), 'ccomplete_outcomeSummary');
    disp(['*** saving: ' fullfile(savepath, ['ccomplete_outcomeSummary_' subjectName '.mat']) ' ***']);
end
    

%% generate averages and export to construct grand averages across mice
cueWindow = [-3 3];
fullWindow = [-6 4];

subSubData = struct(...
    'avg', [],...
    'xData', []...
    );

subData = struct(...
    'licks', subSubData,...
    'photometry', subSubData...
    );

data = struct(...
    'cue', subData,... % cue window
    'full', subData... % full window
    );

% cue
varargin = {'window', cueWindow, 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
avgData = eventAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 'Port1In', varargin{:});
data.cue.licks.avg(1, :, :) = avgData.Avg'; % put trial condition last in matrix
data.cue.licks.xData(1, :, :) = avgData.xData'; 

avgData = phAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 1, 'FluorDataField', 'ZS', 'window', cueWindow, 'zeroTimes', TE.Cue);
data.cue.photometry.avg(1, :, :) = avgData.Avg';
data.cue.photometry.xData(1, :, :) = avgData.xData'; 

% full
varargin = {'window', fullWindow, 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
avgData = eventAverageFromTE(TE, trialsByType, 'Port1In', varargin{:});
data.full.licks.avg(1, :, :) = avgData.Avg'; % put trial condition last in matrix
data.full.licks.xData(1, :, :) = avgData.xData'; 

avgData = phAverageFromTE(TE, trialsByType, 1, 'FluorDataField', 'ZS', 'window', fullWindow, 'zeroTimes', TE.Us);
data.full.photometry.avg(1, :, :) = avgData.Avg';
data.full.photometry.xData(1, :, :) = avgData.xData'; 

ensureFigure('forGrandAverages', 1); 
subplot(2,2,1); plot(squeeze(data.cue.licks.xData(1,:,:)), squeeze(data.cue.licks.avg(1,:,:))); set(gca, 'XLim', cueWindow);
subplot(2,2,3); plot(squeeze(data.cue.photometry.xData(1,:,:)), squeeze(data.cue.photometry.avg(1,:,:))); set(gca, 'XLim', cueWindow);
subplot(2,2,2); plot(squeeze(data.full.licks.xData(1,:,:)), squeeze(data.full.licks.avg(1,:,:))); set(gca, 'XLim', fullWindow);
subplot(2,2,4); plot(squeeze(data.full.photometry.xData(1,:,:)), squeeze(data.full.photometry.avg(1,:,:))); set(gca, 'XLim', fullWindow);
if saveOn
    save(fullfile(savepath, ['averages_' subjectName '.mat']), 'data');
    disp(['*** saving: ' fullfile(savepath, ['averages_' subjectName '.mat']) ' ***']);
end

    
    %% dFF vs reward licks scatter plot
    ensureFigure('phVSLicks_Scatter', 1);
    scatter(TE.csLicks.count(trialsByType{1}) + (rand(length(find(trialsByType{1})), 1) - 0.5) * .1, TE.phPeak_us.data(trialsByType{1}), 'b'); hold on;
    scatter(TE.csLicks.count(trialsByType{4}) + (rand(length(find(trialsByType{4})), 1) - 0.5) * .1, TE.phPeak_us.data(trialsByType{4}), 'r');
        set(gca, 'FontSize', 12); xlabel('cue licks (jittered)'); ylabel('phReward (ZS-avg)');
    set(gca, 'XLim', [0 30]);
if saveOn
    saveas(gcf, fullfile(savepath, 'dFF_vs_licks.fig'));
    saveas(gcf, fullfile(savepath, 'dFF_vs_licks.jpg'));        
end




    %% dFF vs cue licks scatter plot
    ensureFigure('phVS_cueLicks_Scatter_2', 1);
    hrtrials = find(trialsByType{1});
%     hrtrials = hrtrials([1:300 302:end]); % kludge for ChAT_42, NaN in
%  lick count
    lrtrials = find(trialsByType{4});
    scatter(TE.csLicks.count(hrtrials), TE.phPeak_us.data(hrtrials), 'k'); hold on;
    scatter(TE.csLicks.count(trialsByType{4}), TE.phPeak_us.data(lrtrials), 'r');
%     fo = fitoptions('StartPoint', [-5e-4, .005]);
%     fob = fit(TE.csLicks.count(hrtrials), TE.phPeak_us.data(hrtrials), 'poly1', 'options', fo);
%     plot(fob);
    set(gca, 'FontSize', 12); xlabel('cue licks (jittered)'); ylabel('phReward (ZS-avg)');
    set(gca, 'XLim', [0 10]);
% if saveOn
%     saveas(gcf, fullfile(savepath, 'dFF_vs_cueLicks.fig'));
%     saveas(gcf, fullfile(savepath, 'dFF_vs_cueLicks.jpg'));        
% end
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
    


    
