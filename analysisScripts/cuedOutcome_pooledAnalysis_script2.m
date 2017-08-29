% cuedOutcome_Odor_complete summary analysis script 8/28/17
% includes separate analysis for phasic and sustained components to odor
% cue response

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

%% 
savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\Pooled_170828';
ensureDirectory(savepath);
%%
saveOn = 1;
%% 
ensureFigure('CuedOutcome_SummaryData_scatterPlots', 1); 

subplot(2, 2, 1, 'FontSize', 12, 'LineWidth', 1);
scatter(sumData.phCue_phasic_low.avg,sumData.phCue_phasic_high.avg);
addUnityLine;
xlabel('Low Value (Z Score)'); ylabel('High Value (Z Score)'); title('Phasic');

subplot(2, 2, 2, 'FontSize', 12, 'LineWidth', 1);
scatter(sumData.cueLicks_low.avg,sumData.cueLicks_high.avg);
addUnityLine;
xlabel('Low Value (licks/s)'); ylabel('High Value (licks/s)'); title('Anitcipatory Licking');

subplot(2, 2, 3, 'FontSize', 12, 'LineWidth', 1);
scatter(sumData.phCue_sustained_low.avg,sumData.phCue_sustained_high.avg);
addUnityLine;
xlabel('Low Value (Z Score)'); ylabel('High Value (Z Score)'); title('Sustained');

if saveOn
    saveas(gcf, fullfile(savepath, 'CuedOutcome_SummaryData_scatterPlots.fig'));
    saveas(gcf, fullfile(savepath, 'CuedOutcome_SummaryData_scatterPlots.jpg'));
    saveas(gcf, fullfile(savepath, 'CuedOutcome_SummaryData_scatterPlots.epsc'));    
end