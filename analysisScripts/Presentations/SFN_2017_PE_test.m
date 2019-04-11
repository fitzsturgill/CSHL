% SFN 2017 script 
% testing of unsigned prediction error hypothesis

%% desktop
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\Meetings\SFN_2017\PE_Hypothesis';
savepath = 'C:\Users\Fitz\Dropbox\KepecsLab\_Fitz\Meetings\SFN_2017\PE_Hypothesis\';
saveOn = 1;

%% Snippet to make lick histogram for graded value task from  DC_26
figSize = [4.5 3];
% Snippet to make uncued reward and punishment responses from ChAT_26
    load('Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\DC_26\TE.mat');
    cuedOutcome_Conditions;
 % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-6 4], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};

    ensureFigure('Lick_Hist_CCN', 1); axes;

    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 4 7]), 'Port1In', varargin{:});
    set(hl, 'LineWidth', 2);
    legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'Interpreter', 'tex'); legend('boxoff');
    set(gca, 'XLim', [-6 4], 'YLim', [-2 15]);
    addStimulusPatch(gca, [-3 -2 10], '', [0.7 0.7 0.7]) 
    addStimulusPatch(gca, [-0.1 0.1], '');
    ylabel('Licks (Hz)'); xlabel('Time from reward (s)');
    formatFigurePoster(figSize);
    if saveOn    
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN'), 'jpeg');
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN'), 'epsc');
        saveas(gcf, fullfile(savepath, 'Lick_Hist_CCN'), 'fig');        
    end
    
    
%% uncued reward and punishment responses for CBF and dopamine
figSize = [3 2];
% Snippet to make uncued reward and punishment responses from ChAT_26
    load('Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_26\TE.mat');
    cuedOutcome_Conditions;
    saveName = ['Reward_ChAT'];
    ensureFigure(saveName, 1); 
    [ha, hl] = phPlotAverageFromTE(TE, rewardTrials & uncuedTrials, 1,...
        'window', [-1 3], 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'cmap', [171 55 214]/256);
    set(hl, 'LineWidth', 2);
    set(gca, 'XLim', [-1 3], 'YLim', [-1 2], 'YTick', [-1 0 1]);
    addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7]) 
%     ylabel('Fluor. (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from reinforcement (s)');
    formatFigurePoster(figSize);    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    
% Snippet to make uncued reward and punishment responses from ChAT_26
    saveName = ['Punish_ChAT'];
    ensureFigure(saveName, 1); 
    [ha, hl] = phPlotAverageFromTE(TE, punishTrials & uncuedTrials, 1,...
        'window', [-1 3], 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'cmap', [171 55 214]/256);
    set(hl, 'LineWidth', 2);
    set(gca, 'XLim', [-1 3], 'YLim', [-1 2], 'YTick', [-1 0 1]);
    addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7]) 
%     ylabel('Fluor. (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from reinforcement (s)');
    formatFigurePoster(figSize);    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    
    
% Snippet to make uncued reward and punishment responses from DC_26
    load('Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\DC_26\TE.mat');
    cuedOutcome_Conditions;
    saveName = ['Reward_DAT'];
    ensureFigure(saveName, 1); 
    [ha, hl] = phPlotAverageFromTE(TE, rewardTrials & uncuedTrials, 1,...
        'window', [-1 3], 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'cmap', [237 125 49]/256);
    set(hl, 'LineWidth', 2);
    set(gca, 'XLim', [-1 3], 'YLim', [-5 16], 'YTick', [0 10]); %, 'YTickLabel', [0 1]);
    addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7]) 
%     ylabel('Fluor. (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from reinforcement (s)');
    formatFigurePoster(figSize);    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    
    saveName = ['Punish_DAT'];
    ensureFigure(saveName, 1); 
    [ha, hl] = phPlotAverageFromTE(TE, punishTrials & uncuedTrials, 1,...
        'window', [-1 3], 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'cmap', [237 125 49]/256);
    set(hl, 'LineWidth', 2);
    set(gca, 'XLim', [-1 3], 'YLim', [-3 3], 'YTick', [-2 0 2]);%, 'YTickLabel', {'-2', '0', '  2'});
    addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7]) 
%     ylabel('Fluor. (\fontsize{20}\sigma\fontsize{16}-baseline)'); xlabel('Time from reinforcement (s)');
    formatFigurePoster(figSize);    
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    

%% Cs and Us scatter plots
files = {...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_34', 'summary2_ChAT_34.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_35', 'summary2_ChAT_35.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_39', 'summary2_ChAT_39.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_42', 'summary2_ChAT_42.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_26', 'summary2_ChAT_26.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_32', 'summary2_ChAT_32.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_37', 'summary2_ChAT_37.mat'}; % this one last because of weird emergence and dissapearence of reward responses

for counter = 1:size(files, 1)
    fileName = files{counter, 2}(1:end-4);
    temp = load(fullfile(files{counter, 1}, files{counter, 2}));
    if counter == 1
        sumData = temp.cComplete_summary2; % copy it first time, then add
        sumData.filename = {fileName};
    else
        thisData = temp.cComplete_summary2;
        fnames = fieldnames(thisData)';
        sumData.filename{end + 1} = fileName;
        for f = fnames
            sumData.(f{:}).avg(end + 1) = thisData.(f{:}).avg;
            sumData.(f{:}).n(end + 1) = thisData.(f{:}).n;
            sumData.(f{:}).std(end + 1) = thisData.(f{:}).std;
            sumData.(f{:}).sem(end + 1) = thisData.(f{:}).sem;            
        end
    end
end

% Us scatter plot 
figSize = [4,4];
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\Adam\';
saveName = 'CuedOutcome_Us_scatterPlot';
ensureFigure(saveName, 1); 

% {0.9258, 0.4883, 0.1914} % orange
% {0.6680, 0.2148, 0.8359} % purple
% scatter(sumData.phReward_mean_bl.avg,sumData.phPunish_mean_bl.avg, 42, [0.6680, 0.2148, 0.8359], 'filled');
errorbar(sumData.phReward_mean.avg,sumData.phPunish_mean.avg, -sumData.phPunish_mean.sem, sumData.phPunish_mean.sem,...
    -sumData.phReward_mean.sem, sumData.phReward_mean.sem, '.', 'Color', mycolors('chat'));

% xlabel('Reward (\fontsize{20}\sigma\fontsize{16}-baseline)'); ylabel('Punish (\fontsize{20}\sigma\fontsize{16}-baseline)'); 
% xlabel('Reward (\sigma-baseline)'); ylabel('Punish (\sigma-baseline)'); 
setXYsymmetric; addOrginLines(gca, [0 0 0]);
set(gca, 'XTick', [-1 0 1], 'YTick', [-2 0 2]);
formatFigurePoster([3 3]);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    

% cs scatter plot
saveName = 'CuedOutcome_Cs_scatterPlot_Final';
ensureFigure(saveName, 1); 

scatter(sumData.phCue_phasic_low.avg,sumData.phCue_phasic_high.avg, 36, [0.6680, 0.2148, 0.8359], 'filled');
errorbar(sumData.phCue_phasic_low.avg,sumData.phCue_phasic_high.avg, -sumData.phCue_phasic_high.sem, sumData.phCue_phasic_high.sem,...
    -sumData.phCue_phasic_low.sem, sumData.phCue_phasic_low.sem, '.', 'Color', mycolors('chat'));
% xlabel('Low Value (\fontsize{20}\sigma\fontsize{16}-baseline)'); ylabel('High Value (\fontsize{20}\sigma\fontsize{16}-baseline)');

set(gca, 'YLim', [0 1.2]); set(gca, 'XLim', [0 1.2]);
set(gca, 'XTick', [0 0.5 1], 'YTick', [0 0.5 1]);
formatFigurePoster([3 3], '');
% set(gca, 'Position', [0.2 0.2 0.9 0.9]);
addUnityLine(gca, [0 0 0]);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    

%% reward bar graph
files = {...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_34', 'ccomplete_outcomeSummary_ChAT_34.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_35', 'ccomplete_outcomeSummary_ChAT_35.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_39', 'ccomplete_outcomeSummary_ChAT_39.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_42', 'ccomplete_outcomeSummary_ChAT_42.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_26', 'ccomplete_outcomeSummary_ChAT_26.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_32', 'ccomplete_outcomeSummary_ChAT_32.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_37', 'ccomplete_outcomeSummary_ChAT_37.mat'}; % this one last because of weird emergence and dissapearence of reward responses

for counter = 1:size(files, 1)
    fileName = files{counter, 2}(1:end-4);
    temp = load(fullfile(files{counter, 1}, files{counter, 2}));
    if counter == 1
        sumOutcomeData = temp.ccomplete_outcomeSummary; % copy it first time, then add
        sumOutcomeData(1).filename = {fileName};
        fnames = fieldnames(temp.ccomplete_outcomeSummary)';     % field names must all be the same, (avg, n, std, etc.)       
    else
        thisData = temp.ccomplete_outcomeSummary;
        sumOutcomeData(1).filename{end + 1} = fileName;
        for type = 1:9
            for f = fnames
                sumOutcomeData(type).(f{:})(end + 1) = thisData(type).(f{:});
            end
        end          
    end
end
saveName = 'cuedOutcome_reward_barGraph';
ensureFigure(saveName, 1);
offset = 0.1;
bar(1, mean(sumOutcomeData(1).avg), 'b'); hold on; 
bar(2, mean(sumOutcomeData(4).avg), 'r'); 
bar(3, mean(sumOutcomeData(7).avg), 'g'); 
errorbar(1, mean(sumOutcomeData(1).avg), mean(sumOutcomeData(1).avg)/sqrt(length(sumOutcomeData(1).avg)), 'b.', 'linewidth', 2);
errorbar(2, mean(sumOutcomeData(4).avg), mean(sumOutcomeData(4).avg)/sqrt(length(sumOutcomeData(4).avg)), 'r.', 'linewidth', 2);
errorbar(3, mean(sumOutcomeData(7).avg), mean(sumOutcomeData(7).avg)/sqrt(length(sumOutcomeData(7).avg)), 'g.', 'linewidth', 2);
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'High V.', 'Low V.', 'Uncued'}, 'YLim', [0 1.2], 'XLim', [0.5 3.5]);
ylabel('Fluor. (\Delta \sigma-baseline)');
formatFigurePoster([5 4]);   
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    

%% combined cue and licking for ChAT_42

    load('Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_42\TE.mat');
    cuedOutcome_Conditions;
    saveName = 'CuedOutcome_Cue_combined';
    ensureFigure(saveName, 1); 
        varargin = {'window', [-4 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'b', 'r', 'g'}};
    lickAvg = eventAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 'Port1In', varargin{:});
    ax = axes; hold on; yyaxis right
    ll = plot(lickAvg.xData, lickAvg.Avg(1,:), '--k', 'LineWidth', 2);
    plot(lickAvg.xData, lickAvg.Avg(1,:), '--b', lickAvg.xData, lickAvg.Avg(2,:), '--r', lickAvg.xData, lickAvg.Avg(3,:), '--g', 'LineWidth', 2);
    ax.YColor = [0 0 0]; ylabel('Lick rate (Hz)');
    yyaxis left; ax.YColor = [0 0 0]; 
    [ha, hl] = phPlotAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 1,...
        'window', [-4 3], 'linespec', {'b', 'r', 'g'}, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue, 'alpha', 0);
    set(hl, 'LineWidth', 2);    
    set(gca, 'XLim', [-4 3], 'YLim', [-0.2 1.5]);
    addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 1) 
    legend([hl ll(1)], {'\color{blue} high value', '\color{red} low value', '\color{green} uncued', ...
        'licking'}, 'Location', 'northwest', 'FontSize', 20, 'Interpreter', 'tex'); legend('boxoff');
    ylabel('Fluor. (\sigma-baseline)'); xlabel('Time from cue (s)');
    formatFigurePoster([6 4], '', 24);    

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    
    
%%
saveName = 'CuedOutcome_Reward';
ensureFigure(saveName, 1); 
[ha, hl] = phPlotAverageFromTE(TE, {trialsByType{1}, trialsByType{4}, trialsByType{7}}, 1,...
    'window', [-2 2], 'linespec', {'b', 'r', 'g'}, 'FluorDataField', 'ZS', 'alpha', 0);

set(hl, 'LineWidth', 2);
set(gca, 'XLim', [-2 2], 'YLim', [-1.2 2.5]);
addStimulusPatch(gca, [-0.1 0.1], '', [0.8 0.8 0.8], 1);
legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'Location', 'northwest', 'FontSize', 20, 'Interpreter', 'tex'); legend('boxoff');    
ylabel('Fluor. (\sigma-baseline)'); xlabel('Time from reward (s)');
formatFigurePoster([5 4]);   
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    

%% acquisition- from McKnight Poster
% load('C:\Users\Adam\Dropbox\KepecsLab\_Fitz\McKnight\Poster\McKnight.mat');

