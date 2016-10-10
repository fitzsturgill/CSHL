
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
savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete_161007\ChAT_34';
%%
truncateSessionsFromTE(TE, 'init');

%% Make tiled array of licks in receipt of reward to visualize satiation/ lapsing behavior towards end of each session
ensureFigure('RewardLickRate_crossSessions', 1);
for tei = 1:length(TE)
    rewardTrials = filterTE(TE(tei), 'trialOutcome', 1);
    subplot(ceil(sqrt(length(TE))),ceil(sqrt(length(TE))),tei); plot(find(rewardTrials), smooth(TE(tei).usLicks.rate(rewardTrials), 5)); hold on; 
    plot(find(rewardTrials), [0; diff(TE(tei).sessionIndex(rewardTrials))] * 10);
    ylabel('Lick/s, Us'); xlabel('trial #'); textBox(TE(tei).filename{1}(1:7));
    set(gca, 'YLim', [0 10]);
end
saveas(gcf, fullfile(savepath, 'RewardLickRate_crossSessions.fig'));

%% make tiled array of antic. licks for low and high value odors vs trial number
smoothFactor = 11;
ensureFigure('AnticLickRate_crossSessions', 1);
for tei = 1:length(TE)
    highTrials = filterTE(TE(tei), 'trialType', 1:3, 'reject', 0);
    lowTrials = filterTE(TE(tei), 'trialType', 4:6, 'reject', 0);
    rewardTrials = filterTE(TE(tei), 'trialOutcome', 1);    
    subplot(ceil(sqrt(length(TE))), ceil(sqrt(length(TE))), tei); plot(find(highTrials), smooth(TE(tei).csLicks.rate(highTrials), smoothFactor), 'b.'); hold on; 
    plot(find(lowTrials), smooth(TE(tei).csLicks.rate(lowTrials), smoothFactor), 'r.');
    plot(find(rewardTrials), smooth(TE(tei).usLicks.rate(rewardTrials), smoothFactor), 'k.')    
    plot(1:length(TE(tei).trialNumber), [0; diff(TE(tei).sessionIndex)] * 10, 'g');
    ylabel('Lick/s, Us'); xlabel('trial #'); textBox(TE(tei).filename{1}(1:7));
    set(gca, 'YLim', [0 10]);
end
saveas(gcf, fullfile(savepath, 'AnticLickRate_crossSessions.fig'));

%% figure out at which trial mouse lapses/ becomes sated, go through all 6 mice and all the sessions
% to find out num of trials in session: max(find(filterTE(TE(mi), 'sessionIndex', si))) - min(find(filterTE(TE(mi), 'sessionIndex', si)));
si = 1; % session index [1:6 8:10]
mi = 1; % mouse index
trials = filterTE(TE(mi), 'sessionIndex', si, 'trialType', 1);
inputArgs = {'startField', 'PreCsRecording',...
    'zeroField', 'Us',...
    'endField', 'PostUsRecording',...
    'trialNumbering', 'singleSession'};
ensureFigure('lickRaster', 1); axes; eventRasterFromTE(TE(mi), trials, 'Port1In', inputArgs{:});
%% determined based upon inspection of when antic. licking for high value cue falls off, see above
TrialCutoffs{1} = [139 160 230 160 170 190 254 165]; % ChAT_39
% TrialCutoffs{1} = [110 125 80 110 125 70 200 150];
% TrialCutoffs{2} = [87 140 100 85 115 115 120];
% TrialCutoffs{3} = [130 115 110 150 202 125];
% TrialCutoffs{4} = [75 85 90 100 85 130];
% TrialCutoffs{5} = [100 170 100 115 160 200 204];
% TrialCutoffs{6} = [102 110 85 90 140 100 100 80 95 135]; % ChAT_26
%% use TrialCutoffs to add reject field to TE
TE = ensureField(TE, 'reject', 'mat');
for counter = 1:length(TrialCutoffs)
    if ~isempty(TrialCutoffs{counter})
        nSessions = length(unique(TE(counter).sessionIndex));
        if nSessions == length(TrialCutoffs{counter})
            reject = ones(size(TE(counter).trialNumber));
            for i = 1:nSessions
                trialNumbers = [];
                trialNumbers(:,1) = TE(counter).trialNumber;
                include = (TE(counter).sessionIndex == i) & (trialNumbers <= TrialCutoffs{counter}(i));
                reject(include) = 0;
            end
            TE(counter).reject = reject;
        end
    end
end

%% generate trial lookups for different combinations of conditions
    tei = 1;
    includedSessions = [1:8];
    highValueTrials = filterTE(TE(tei), 'trialType', 1:3, 'reject', 0, 'sessionIndex', includedSessions);
    lowValueTrials = filterTE(TE(tei), 'trialType', 4:6, 'reject', 0, 'sessionIndex', includedSessions);
    rewardTrials = filterTE(TE(tei), 'trialOutcome', 1, 'reject', 0, 'sessionIndex', includedSessions);
    punishTrials = filterTE(TE(tei), 'trialOutcome', 2, 'reject', 0, 'sessionIndex', includedSessions);
    omitTrials = filterTE(TE(tei), 'trialOutcome', 3, 'reject', 0, 'sessionIndex', includedSessions);
    trialTypes = 1:9;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE(tei), 'trialType', trialTypes(counter), 'reject', 0, 'sessionIndex', includedSessions);
    end
    %% plot photometry averages
    ensureFigure('Photometry_Averages', 1); 
    pm = [3 2];
    
    % - 6 0 4
    subplot(pm(1), pm(2), 1); [ha, hl] = phPlotAverageFromTE(TE(end), trialsByType([1 3 7]), 1); %high value, reward
    legend(hl, {'hival, rew', 'hival, omit', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('high value'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 2); [ha, hl] = phPlotAverageFromTE(TE(end), trialsByType([5 6 8]), 1); % low value, punish
    legend(hl, {'loval, pun', 'loval, omit', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('low value'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 3); [ha, hl] = phPlotAverageFromTE(TE(end), trialsByType([1 4 7]), 1); % reward, varying degrees of expectation
    legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('reward all'); ylabel('dF/F'); xlabel('time from reinforcement (s)');     

    subplot(pm(1), pm(2), 4); [ha, hl] = phPlotAverageFromTE(TE(end), trialsByType([5 2 8]), 1); % punishment, varying degrees of expectation
    legend(hl, {'loval, pun', 'hival, pun', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('punish all'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 5); [ha, hla] = phPlotAverageFromTE(TE(end), {lowValueTrials, highValueTrials}, 1,...
        'window', [-6 0], 'linespec', {'m', 'g'}); hold on;
    
    subplot(pm(1), pm(2), 5); [ha, hl] = phPlotAverageFromTE(TE(end), {rewardTrials, punishTrials, omitTrials}, 1,...
        'window', [0 4], 'linespec', {'b', 'r', 'k'});
    hl = [hla hl];
    legend(hl, {'loval', 'hival', 'rew', 'pun', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Balazs'); ylabel('dF/F'); xlabel('time from reinforcement (s)'); 
    
    saveas(gcf, fullfile(savepath, 'phAverages_CuedOutcome_Odor_Complete_ChAT_26.fig'));   
    
    %% plot photometry rasters
    ensureFigure('phRasters', 1); 
    subplot(1,2,1); imshow(TE(end).Photometry.data(1).dFF(trialsByType{1}, :), [-.01 .01]); 
    colormap default; title('hival, reward'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2); imshow(TE(end).Photometry.data(1).dFF(trialsByType{5}, :), [-.01 .01]); 
    colormap default; title('loval, punish'); xlabel('time from reinforcement (s)');
     saveas(gcf, fullfile(savepath, 'phRasters.fig'));  
    
    %% plot photometry rasters- test phRasterFromTE
    ensureFigure('phRastersFromTE', 1); 
    subplot(1,2,1); phRasterFromTE(TE(end), trialsByType{1}, 1);
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2); phRasterFromTE(TE(end), trialsByType{5}, 1);
    title('loval, punish'); xlabel('time from reinforcement (s)');     

    %% off conditions
    %% plot photometry rasters
    ensureFigure('phRasters_off', 1); 
    subplot(1,2,1); imshow(TE(end).Photometry.data(1).dFF(trialsByType{2}, :), [-.01 .01]); 
    colormap default; title('hival, punish'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2); imshow(TE(end).Photometry.data(1).dFF(trialsByType{4}, :), [-.01 .01]); 
    colormap default; title('loval, reward'); xlabel('time from reinforcement (s)');     

    %% plot photometry rasters 2 lab meeting
    ensureFigure('phRasters_hival', 1); 
    prt = trialsByType{1};
    prcd = TE(end).Photometry.data(1).dFF(prt, :);
    subplot(1,2,1); 
    image('Xdata', [-6 4], 'YData', [1 find(length(prt))],...
        'CData', prcd, 'CDataMapping', 'Scaled', 'Parent', gca);
    set(gca, 'CLim', [-0.01 .01], 'YDir', 'Reverse');
    set(gca, 'XLim', [-6 4]);
    set(gca, 'FontSize', 14)
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2);
    eventRasterFromTE(TE(end), trialsByType{1}, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]); 
    set(gca, 'FontSize', 14)
    saveas(gcf, fullfile(savepath, 'phRasters_hiVal.fig'));  
    
    %% plot photometry rasters lab meeting low value
    ensureFigure('phRasters_lowVal', 1); 
    prt = trialsByType{5};
    prcd = TE(end).Photometry.data(1).dFF(prt, :);
    subplot(1,2,1); 
    image('Xdata', [-6 4], 'YData', [1 find(length(prt))],...
        'CData', prcd, 'CDataMapping', 'Scaled', 'Parent', gca);
    set(gca, 'CLim', [-0.01 .01], 'YDir', 'Reverse');
    set(gca, 'XLim', [-6 4]);
    set(gca, 'FontSize', 14)
    title('loVal, punish'); xlabel('time from reinforcement (s)'); 
    subplot(1,2,2);
    eventRasterFromTE(TE(end), prt, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]);        
    set(gca, 'FontSize', 14)
    saveas(gcf, fullfile(savepath, 'phRasters_lowVal.fig'));  
    %% plot lick rasters
    ensureFigure('lickRasters', 1);
    subplot(1,2,1);
    eventRasterFromTE(TE(end), trialsByType{1}, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('hival, reward'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]);
    
    subplot(1,2,2);
    eventRasterFromTE(TE(end), trialsByType{5}, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('lowval, punish'); xlabel('time from reinforcement (s)'); 
    set(gca, 'XLim', [-6 4]);
    
    
