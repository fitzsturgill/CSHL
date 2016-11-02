%%
files = {...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete_clean\ChAT_34', 'summary_ChAT_34.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete_clean\ChAT_35', 'summary_ChAT_35.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_39', 'summary_ChAT_39.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_42', 'summary_ChAT_42.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_26', 'summary_ChAT_26.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_32', 'summary_ChAT_32.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_37', 'summary_ChAT_37.mat'}; % this one last because of weird emergence and dissapearence of reward responses

% concatenate summary data
allData=struct();
for counter = 1:size(files, 1)
    fileName = files{counter, 2}(1:end-4);
    temp = load(fullfile(files{counter, 1}, files{counter, 2}));
    sumData = temp.cComplete_summary;
    fnames = fieldnames(sumData)';
    allData.filename{counter, 1} = fileName;
    for f = fnames
        allData.(f{:})(counter, :) = sumData.(f{:});
    end
end

%%
%          filename: {5x1 cell}
%          phCue_CV: [5x9 double]
%         phCue_avg: [5x9 double]
%      phOutcome_CV: [5x9 double]
%     phOutcome_avg: [5x9 double]
%      cueLicks_low: [5x1 double]
%     cueLicks_high: [5x1 double]
%       rewardLicks: [5x1 double]

savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\Pooled_161102';
ensureDirectory(savepath);

%%
ensureFigure('cuedOutcome_pooledAnalysis', 1); 
subplot(2,2,1);
scatter(allData.phOutcome_avg(:,7), allData.phOutcome_avg(:,8));
addUnityLine;
xlabel('reward dFF'); ylabel('punish dFF');

subplot(2,2,2);
scatter(allData.phOutcome_CV(:,7), allData.phOutcome_CV(:,8));
addUnityLine;
xlabel('reward dFF ~Zscored'); ylabel('punish dFF ~Zscored');

subplot(2,2,3);
scatter(allData.cueLicks_low ./ allData.rewardLicks, allData.cueLicks_high ./ allData.rewardLicks);
addUnityLine;
xlabel('low:reward lick ratio'); ylabel('high:reward lick ratio');

subplot(2,2,4);
scatter(allData.phOutcome_avg(:,4), allData.phOutcome_avg(:,1));
addUnityLine;
xlabel('hival reward dFF'); ylabel('loval reward dFF');

saveas(gcf, fullfile(savepath, 'pooledAnalysis_scatterPlots.fig'));
saveas(gcf, fullfile(savepath, 'pooledAnalysis_scatterPlots.jpg'));
        