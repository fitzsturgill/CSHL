
saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_CuedOutcome_Odor_Complete(sessions);
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial');

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
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
% basepath = uigetdir;
basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\initial_sessions';
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);
disp(['*** ensuring directory: ' savepath ' ***']);
%% truncate sessions according to when mouse becomes sated
rewardTrialsTrunc = filterTE(TE, 'trialOutcome', 1);
truncateSessionsFromTE(TE, 'init', 'usLicks', rewardTrialsTrunc);
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
    %% plot photometry averages
    h=ensureFigure('Photometry_Averages', 1); 
%     mcLandscapeFigSetup(h);

    pm = [1 1];
    
    % - 6 0 4
%     subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 3 7]), 1, 'FluorDataField', 'ZS'); %high value, reward
%     legend(hl, {'hival, rew', 'hival, omit', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
%     title('high value'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); textBox(TE.filename{1}(1:7));


    subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); [ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
        'window', [-6 0], 'linespec', {'m', 'g'}, 'FluorDataField', 'ZS'); hold on;
    
    subplot(pm(1), pm(2), 1); [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, omitTrials}, 1,...
        'window', [0 4], 'linespec', {'b', 'r', 'k'}, 'FluorDataField', 'ZS');
    hl = [hla hl];
    legend(hl, {'loval', 'hival', 'rew', 'pun', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
     ylabel('Z Score'); xlabel('time from reinforcement (s)'); 
    

    
if saveOn    
    saveas(gcf, fullfile(savepath, 'phAverages.fig'));
    saveas(gcf, fullfile(savepath, 'phAverages.jpg'));
end




  



%% plot photometry rasters, high value
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
    