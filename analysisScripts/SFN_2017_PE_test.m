% SFN 2017 script 
% testing of unsigned prediction error hypothesis

%% desktop
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SFN_2017\PE_Hypothesis';
saveOn = 1;

%% Snippet to make lick histogram for graded value task from or DC_26
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
    saveas(gcf, fullfile(savepath, [saveNamse '.epsc']));   
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
scatter(sumData.phReward_mean_bl.avg,sumData.phPunish_mean_bl.avg, 42, [0.6680, 0.2148, 0.8359], 'filled');
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

%%
open('C:\Users\Adam\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\Adam\CuedOutcome_Cue_combined.fig');
% formatFigurePoster([5 4]);
