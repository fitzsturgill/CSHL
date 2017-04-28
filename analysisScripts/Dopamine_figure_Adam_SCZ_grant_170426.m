
saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_CuedOutcome_Odor_Complete(sessions);


% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRGECO1a (ch2)
channels=[]; dFFMode = {};
if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'expFit';
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple'; % or 'simple'
end

TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'channels', channels, 'downsample', 305); % 61 for ph raster


%% extract peak trial dFF responses to cues and reinforcement and lick counts

for channel = channels
    TE.phPeak_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 2], TE.Cue, 'method', 'mean');
    TE.phPeak_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 1], TE.Us, 'method', 'mean');

    TE.csLicks = countEventFromTE(TE, 'Port1In', [-2 0], TE.Us);
    TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.Us);
end
%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = uigetdir;
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
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




%% generate trial lookups for different combinations of conditions
    validTrials = filterTE(TE, 'reject', 0);
    highValueTrials = filterTE(TE, 'trialType', 1:3, 'reject', 0);
    lowValueTrials = filterTE(TE, 'trialType', 4:6, 'reject', 0);
    uncuedTrials = filterTE(TE, 'trialType', 7:9, 'reject', 0);    
    rewardTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    punishTrials = filterTE(TE, 'trialOutcome', 2, 'reject', 0);    
    omitTrials = filterTE(TE, 'trialOutcome', 3, 'reject', 0);
    phRasterTrials = filterTE(TE, 'sessionIndex', [1 2], 'reject', 0);
    trialTypes = 1:9;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end
    %% plot photometry averages
    for channel = channels
        figname = ['Photometry_Averages_ch' num2str(channel)];
        h=ensureFigure(figname, 1); 
        mcLandscapeFigSetup(h);

        pm = [3 2];

        % - 6 0 4
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 3 7]), channel); %high value, reward
        legend(hl, {'hival, rew', 'hival, omit', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('high value'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); textBox(TE.filename{1}(1:7));

        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 6 8]), channel); % low value, punish
        legend(hl, {'loval, pun', 'loval, omit', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('low value'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 
%% reward all
        channel = 2;
        figname = ['Photometry_Averages_rewardAll_ch' num2str(channel)];
        h=ensureFigure(figname, 1); 
        subplot(1, 1, 1, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), channel, 'linespec', {'b', 'r', 'k'}, 'window', [-2 4]); % reward, varying degrees of expectation
        legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        ylabel('dF/F'); xlabel('Time from reinforcement (s)');  
        set(gca, 'XLim', [-2 4]);
        formatFigureGRC([4,3]);
            saveas(gcf, fullfile(savepath, [figname '.fig']));
            saveas(gcf, fullfile(savepath, [figname '.epsc']));
            
%%            
        subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 2 8]), channel); % punishment, varying degrees of expectation
        legend(hl, {'loval, pun', 'hival, pun', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('punish all'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 

        %% Balazs
        
        channel = 2;
        figname = ['Photometry_Averages_Balazs_ch' num2str(channel)];
        h=ensureFigure(figname, 1);         
        subplot(1,1, 1, 'FontSize', 12, 'LineWidth', 1); [ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, channel,...
            'window', [-5 0], 'linespec', {'b', 'r'}); hold on;

        subplot(1,1, 1); [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, omitTrials}, channel,...
            'window', [0 4], 'linespec', {'c', 'm', 'g'});
        hl = [hla hl];
        legend(hl, {'loval', 'hival', 'rew', 'pun', 'omit'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        ylabel('dF/F'); xlabel('Time from reinforcement (s)'); 
        formatFigureGRC([4,3]);
%         if saveOn    
            saveas(gcf, fullfile(savepath, [figname '.fig']));
            saveas(gcf, fullfile(savepath, [figname '.epsc']));
%         end
    end

%% graph to pick phasic and sustained analysis windows
ensureFigure('windowPick', 1);
axes;
phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
        'window', [-4 0], 'linespec', {'m', 'g'}); 
    

%% plot lick averages
    h = ensureFigure('Lick_Averages_SCZ', 1);
%     mcPortraitFigSetup(h);
%     mcLandscapeFigSetup(h);
    pm = [1 1];
    
    % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-5 4], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'linespec', {'b', 'r', 'k'}};
    axh = [];

    
    % reward
    axh(end + 1) = subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 4 7]), 'Port1In', varargin{:});
    legend(hl, {'hival', 'loval', 'uncued'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    ylabel('Licks (1/s)'); xlabel('Time from reinforcement (s)');
    
%     % punish
%     axh(end + 1) = subplot(pm(1), pm(2), 3); [ha, hl] = plotEventAverageFromTE(TE, trialsByType([2 5 8]), 'Port1In', varargin{:});
%     legend(hl, {'hival', 'loval', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
%     title('Punish'); ylabel('licks (s)'); xlabel('time r (s)');
%     
%     % neutral
%     axh(end + 1) = subplot(pm(1), pm(2), 4); [ha, hl] = plotEventAverageFromTE(TE, trialsByType([3 6 9]), 'Port1In', varargin{:});
%     legend(hl, {'hival', 'loval', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
%     title('Neutral'); ylabel('licks (s)'); xlabel('time r (s)');
    
%     sameYScale(axh) % match y scaling
formatFigureGRC([4,3]);
if saveOn    
    saveas(gcf, fullfile(savepath, 'lickAverages_SCZ.fig'));
    saveas(gcf, fullfile(savepath, 'lickAverages_SCZ.epsc'));    
end
  

    %% plot photometry rasters
        channel = 2;
        CLimFactor = 3;
        savename = ['phRasters_SCZ_high_vs_low' num2str(channel)];
        h=ensureFigure(savename, 1);
        mcPortraitFigSetup(h);
        set(gca, 'FontSize', 12)
        subplot(1,2,1); phRasterFromTE(TE, trialsByType{1} & phRasterTrials, channel, 'CLimFactor', CLimFactor);
        set(gca, 'FontSize', 12)
        subplot(1,2,2); phRasterFromTE(TE, trialsByType{4} & phRasterTrials, channel, 'CLimFactor', CLimFactor);
        set(gca, 'FontSize', 12)
%         subplot(1,3,3); phRasterFromTE(TE, trialsByType{7} & phRasterTrials, channel, 'CLimFactor', CLimFactor);
            saveas(gcf, fullfile(savepath, [savename '.fig']));
            saveas(gcf, fullfile(savepath, [savename '.epsc']));
    
    
    
    %% 
    for channel = channels
        CLimFactor = 2;
        savename = ['phRastersFromTE_reward_ch' num2str(channel)];
        h=ensureFigure(savename, 1);
        mcPortraitFigSetup(h);
        subplot(1,4,1); phRasterFromTE(TE, trialsByType{1}, channel, 'CLimFactor', CLimFactor);
        title([TE.filename{1}(1:7) ': hival, reward'], 'Interpreter', 'none'); 
        subplot(1,4,2); phRasterFromTE(TE, trialsByType{4}, channel, 'CLimFactor', CLimFactor);
        title('loval, reward'); xlabel('time from reinforcement (s)'); 
        subplot(1,4,3); phRasterFromTE(TE, trialsByType{7}, channel, 'CLimFactor', CLimFactor);
        title('uncued, reward');
        subplot(1,4,4); phRasterFromTE(TE, trialsByType{3}, channel, 'CLimFactor', CLimFactor);
        title('hival, neutral'); xlabel('time from tone (s)'); 
        
        
        
        
        if saveOn
            saveas(gcf, fullfile(savepath, [figname '.fig']));
            saveas(gcf, fullfile(savepath, [figname '.jpg']));  
        end

        savename = ['phRastersFromTE_punish_ch' num2str(channel)];
        h=ensureFigure(savename, 1);
        mcPortraitFigSetup(h);
        subplot(1,4,1); phRasterFromTE(TE, trialsByType{2}, channel, 'CLimFactor', CLimFactor);
        title([TE.filename{1}(1:7) ': hival, punish'], 'Interpreter', 'none'); 
        subplot(1,4,2); phRasterFromTE(TE, trialsByType{5}, channel, 'CLimFactor', CLimFactor);
        title('loval, punish'); xlabel('time from reinforcement (s)'); 
        subplot(1,4,3); phRasterFromTE(TE, trialsByType{8}, channel, 'CLimFactor', CLimFactor);
        title('uncued, punish');
        subplot(1,4,4); phRasterFromTE(TE, trialsByType{6}, channel, 'CLimFactor', CLimFactor);
        title('loval, neutral'); xlabel('time from tone (s)'); 
        if saveOn
            saveas(gcf, fullfile(savepath, [figname '.fig']));
            saveas(gcf, fullfile(savepath, [figname '.jpg']));
        end
    end


%%
    %% plot photometry rasters reward- alternate for lab meeting
    for channel = channels
        CLimFactor = 2;
        savename = ['phRastersAlternate_reward_ch' num2str(channel)];
        h=ensureFigure(savename, 1);
        mcPortraitFigSetup(h);


    %     prcd = TE.Photometry.data(1).dFF(prt, :);
        subplot(1,4,1); 
        eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', 'consecutive',...
            'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
        title('hival, reward'); xlabel('time from reinforcement (s)'); ylabel('trial number');
        set(gca, 'XLim', [-6 4]); 
        set(gca, 'YLim', [0 length(find(trialsByType{1}))]);
        set(gca, 'FontSize', 14)


        subplot(1,4,2); phRasterFromTE(TE, trialsByType{1}, channel, 'CLimFactor', CLimFactor);
        title('hival, reward', 'Interpreter', 'none'); 
            set(gca, 'FontSize', 14)
        subplot(1,4,3); phRasterFromTE(TE, trialsByType{4}, channel, 'CLimFactor', CLimFactor);
        title('loval, reward'); 
            set(gca, 'FontSize', 14)
        subplot(1,4,4); phRasterFromTE(TE, trialsByType{7}, channel, 'CLimFactor', CLimFactor);
        title('uncued, reward');
            set(gca, 'FontSize', 14)


        saveas(gcf, fullfile(savepath, [savename '.fig']));
        saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    end

% %% plot photometry rasters 2 lab meeting
%     h=ensureFigure('phRasters_hival', 1);
%     mcPortraitFigSetup(h);
%     prt = trialsByType{1};
% %     prcd = TE.Photometry.data(1).dFF(prt, :);
%     subplot(1,2,1); 
%     phRasterFromTE(TE, trialsByType{1}, 1, 'CLimFactor', CLimFactor);
% %     image('Xdata', [-6 4], 'YData', [1 find(length(prt))],...
% %         'CData', prcd, 'CDataMapping', 'Scaled', 'Parent', gca);
% %     set(gca, 'CLim', [-0.01 .01], 'YDir', 'Reverse');
%     set(gca, 'XLim', [-6 4]);
%     set(gca, 'YLim', [0 length(find(trialsByType{1}))]);
%     set(gca, 'FontSize', 14)
%     title('hival, reward'); xlabel('time from reinforcement (s)'); 
%     subplot(1,2,2);
%     eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', 'consecutive',...
%         'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
%     title('hival, reward'); xlabel('time from reinforcement (s)'); 
%     set(gca, 'XLim', [-6 4]); 
%     set(gca, 'YLim', [0 length(find(trialsByType{1}))]);
%     set(gca, 'FontSize', 14)
% if saveOn
%     saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_reward.fig'));
%     saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_reward.jpg'));    
% end
%     %% plot photometry rasters lab meeting low value
%     h=ensureFigure('phRasters_lowVal', 1); 
%     mcPortraitFigSetup(h);    
%     prt = trialsByType{5};
% %     prcd = TE.Photometry.data(1).dFF(prt, :);
%     subplot(1,2,1); 
%     phRasterFromTE(TE, trialsByType{5}, 1, 'CLimFactor', CLimFactor);    
% %     image('Xdata', [-6 4], 'YData', [1 find(length(prt))],...
% %         'CData', prcd, 'CDataMapping', 'Scaled', 'Parent', gca);
%     set(gca, 'CLim', [-0.01 .01], 'YDir', 'Reverse');
%     set(gca, 'YLim', [0 length(find(trialsByType{5}))]);   
%     set(gca, 'XLim', [-6 4]);
%     set(gca, 'FontSize', 14)
%     title('loVal, punish'); xlabel('time from reinforcement (s)'); 
%     subplot(1,2,2);
%     eventRasterFromTE(TE, prt, 'Port1In', 'trialNumbering', 'consecutive',...
%         'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
%     title('loVal, punish'); xlabel('time from reinforcement (s)'); 
%     set(gca, 'XLim', [-6 4]);     
%     set(gca, 'YLim', [0 length(find(trialsByType{5}))]);        
%     set(gca, 'FontSize', 14)
% if saveOn
%     saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_punish.fig'));
%     saveas(gcf, fullfile(savepath, 'licks_ph_comp_raster_punish.jpg'));    
% end
%     %% summary statistics
%     cComplete_summary = struct(...
%         'phCue_CV', zeros(1,9),...
%         'phCue_avg', zeros(1,9),...
%         'phOutcome_CV', zeros(1,9),...
%         'phOutcome_avg', zeros(1,9),...
%         'cueLicks_low', 0,...
%         'cueLicks_high', 0,...
%         'rewardLicks', 0);
%     for counter = 1:length(trialTypes)
%         trials = trialsByType{counter};
%         peaks = TE.phPeak_cs.data(trials);
%         cComplete_summary.phCue_CV(counter) = nanmean(peaks) / std(peaks, 'omitnan');
%         cComplete_summary.phCue_avg(counter) = nanmean(peaks);
%         peaks = TE.phPeak_us.data(trials);
%         cComplete_summary.phOutcome_CV(counter) = nanmean(peaks) / std(peaks, 'omitnan');
%         cComplete_summary.phOutcome_avg(counter) = nanmean(peaks); 
%     end
%     cComplete_summary.cueLicks_low = nanmean(TE.csLicks.rate(lowValueTrials));
%     cComplete_summary.cueLicks_high = nanmean(TE.csLicks.rate(highValueTrials));        
%     cComplete_summary.rewardLicks = nanmean(TE.usLicks.rate(rewardTrials));    
% if saveOn
%     save(fullfile(savepath, ['summary_' subjectName '.mat']), 'cComplete_summary');
%     disp(['*** saving: ' fullfile(savepath, ['summary_' subjectName '.mat']) ' ***']);
% end

    
    

    
%     %% dFF vs reward licks scatter plot
%     ensureFigure('phVSLicks_Scatter', 1);
%     scatter(TE.csLicks.count(trialsByType{1}) + (rand(length(find(trialsByType{1})), 1) - 0.5) * .1, TE.phPeak_us.data(trialsByType{1}), 'b'); hold on;
%     scatter(TE.csLicks.count(trialsByType{4}) + (rand(length(find(trialsByType{4})), 1) - 0.5) * .1, TE.phPeak_us.data(trialsByType{4}), 'r');
%         set(gca, 'FontSize', 12); xlabel('cue licks (jittered)'); ylabel('phReward (dFF-avg)');
%     set(gca, 'XLim', [0 30]);
% if saveOn
%     saveas(gcf, fullfile(savepath, 'dFF_vs_licks.fig'));
%     saveas(gcf, fullfile(savepath, 'dFF_vs_licks.jpg'));        
% end




%     %% dFF vs cue licks scatter plot
%     ensureFigure('phVS_cueLicks_Scatter_2', 1);
%     hrtrials = find(trialsByType{1});
% %     hrtrials = hrtrials([1:300 302:end]); % kludge for ChAT_42, NaN in
% %  lick count
%     lrtrials = find(trialsByType{4});
%     scatter(TE.csLicks.count(hrtrials), TE.phPeak_us.data(hrtrials), 'k'); hold on;
%     scatter(TE.csLicks.count(trialsByType{4}), TE.phPeak_us.data(lrtrials), 'r');
% %     fo = fitoptions('StartPoint', [-5e-4, .005]);
% %     fob = fit(TE.csLicks.count(hrtrials), TE.phPeak_us.data(hrtrials), 'poly1', 'options', fo);
% %     plot(fob);
%     set(gca, 'FontSize', 12); xlabel('cue licks (jittered)'); ylabel('phReward (dFF-avg)');
%     set(gca, 'XLim', [0 10]);
% % if saveOn
% %     saveas(gcf, fullfile(savepath, 'dFF_vs_cueLicks.fig'));
% %     saveas(gcf, fullfile(savepath, 'dFF_vs_cueLicks.jpg'));        
% % end
%     %% us vs cs dFF for highValue reward condition
%     ensureFigure('phUSvsCS_Scatter', 1);
%     scatter(TE.phPeak_cs.data(trialsByType{1}), TE.phPeak_us.data(trialsByType{1}), 'b'); hold on;
% if saveOn
%     saveas(gcf, fullfile(savepath, 'us_vs_cs_dFF.fig'));
%     saveas(gcf, fullfile(savepath, 'us_vs_cs_dFF.jpg'));    
% end
%     %% is the trough related to number of licks? 
%     TE.phTrough_us = bpCalcPeak_dFF(TE.Photometry, 1, [0 2], TE.Us, 'method', 'min');
%     ensureFigure('phDip_Scatter', 1);
%     scatter(TE.usLicks.count(rewardTrials) + rand(length(find(rewardTrials)), 1) - 0.5, TE.phTrough_us.data(rewardTrials), 'b'); 
%     
    %%  reversal experiments
    
    %% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks    blocks 2 and 3
    validTrials = filterTE(TE, 'reject', 0);
    Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
    Odor2Trials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);    
    
    rewardTrials = filterTE(TE, 'trialType', 1, 'reject', 0);
    hitTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    missTrials = filterTE(TE, 'trialOutcome', -1, 'reject', 0);
    punishTrials = filterTE(TE, 'trialType', 3, 'reject', 0);    
    neutralTrials = filterTE(TE, 'trialType', [2 4], 'reject', 0);
    block2Trials = filterTE(TE, 'BlockNumber', 2, 'reject', 0);
    block3Trials = filterTE(TE, 'BlockNumber', 3, 'reject', 0);
    csPlusTrials = filterTE(TE, 'trialType', [1 2], 'reject', 0);
    csMinusTrials = filterTE(TE, 'trialType', [3 4], 'reject', 0);    
    trialTypes = 1:4;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end

trialCount = [1:length(TE.filename)]';

%%

%% PH Rasters, Odor 1, Odor 2
CLimFactor = 2;

reversals = find(TE.BlockChange);
for channel = channels

    saveName = [subjectName '_comboRasters_SCZ_ch' num2str(channel)];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);
    

%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,4,1); 
    eventRasterFromTE(TE, Odor1Trials, 'Port1In', 'trialNumbering', 'global',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('Odor A'); ylabel('trial number');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    

    
    subplot(1,4,2); phRasterFromTE(TE, Odor1Trials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
        set(gca, 'FontSize', 14)
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines        
        
    subplot(1,4,3); 
    eventRasterFromTE(TE, Odor2Trials, 'Port1In', 'trialNumbering', 'global',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    title('Odor B');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    
    subplot(1,4,4); phRasterFromTE(TE, Odor2Trials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    set(gca, 'FontSize', 14)
    xlabel('time from cue (s)');         


    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));
    end
end

%% hit vs miss phAverages
    saveName = 'hit_vs_miss_phAvgs_ch2';  
    h=ensureFigure(saveName, 1); 
    
    
    subplot(1, 1, 1, 'FontSize', 12, 'LineWidth', 1); 
    [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 2,...
        'FluorDataField', 'ZS', 'window', [-4, 7], 'linespec', {'c', 'm'}); %high value, reward
    legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    xlabel('time from cue (s)');     set(gca, 'XLim', [-4, 7]);

    formatFigureGRC;
