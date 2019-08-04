% cuedOutcome_pooled_figure

% goals: summary data for cue and outcome responses


DB = dbLoadExperiment('cuedOutcome');

photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
ch = 1;

savePath = fullfile(DB.path, 'pooled');
ensureDirectory(savePath);
figPath = fullfile(DB.path, 'figure');
ensureDirectory(savePath);

% compile cue and outcome responses for high value, low value, and null,
% reward and punish

nAnimals = length(DB.animals);
s2 = struct(...
    'data', [],...
    'avg', zeros(nAnimals, 1),...
    'SEM', zeros(nAnimals, 1)...
    );
cs_pooled = struct(...
    'high', s2,...
    'low', s2,...
    'uncued', s2...
    );
us_pooled = struct(...
    'reward', s2,... % uncued
    'punish', s2...  % uncued
    );

for counter = 1:nAnimals
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);        
    
    thisData = TE.phPeak_cs.data(highValueTrials);
    cs_pooled.high.data{counter} = thisData;
    cs_pooled.high.avg(counter) = nanmean(thisData);
    cs_pooled.high.SEM(counter) = nanSEM(thisData);
    
    thisData = TE.phPeak_cs.data(lowValueTrials);
    cs_pooled.low.data{counter} = thisData;
    cs_pooled.low.avg(counter) = nanmean(thisData);
    cs_pooled.low.SEM(counter) = nanSEM(thisData);
    
    thisData = TE.phPeak_us.data(uncuedTrials & rewardTrials);
    us_pooled.reward.data{counter} = thisData;
    us_pooled.reward.avg(counter) = nanmean(thisData);
    us_pooled.reward.SEM(counter) = nanSEM(thisData);
    
    thisData = TE.phPeak_us.data(uncuedTrials & rewardTrials);
    us_pooled.reward.data{counter} = thisData;
    us_pooled.reward.avg(counter) = nanmean(thisData);
    us_pooled.reward.SEM(counter) = nanSEM(thisData);    
    
    thisData = TE.phPeak_us.data(uncuedTrials & punishTrials);
    us_pooled.punish.data{counter} = thisData;
    us_pooled.punish.avg(counter) = nanmean(thisData);
    us_pooled.punish.SEM(counter) = nanSEM(thisData);        
end

%%
save(fullfile(savePath, 'us_pooled.mat'), 'us_pooled');
disp(['*** saving: ' fullfile(savePath, 'us_pooled.mat') ' ***']);
save(fullfile(savePath, 'cs_pooled.mat'), 'cs_pooled');
disp(['*** saving: ' fullfile(savePath, 'cs_pooled.mat') ' ***']);

%% us scatter plot

% Us scatter plot 
figSize = [1 1];
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\Adam\';
saveName = 'CuedOutcome_Us_scatterPlot';
ensureFigure(saveName, 1); 

% {0.9258, 0.4883, 0.1914} % orange
% {0.6680, 0.2148, 0.8359} % purple
% scatter(sumData.phReward_mean_bl.avg,sumData.phPunish_mean_bl.avg, 42, [0.6680, 0.2148, 0.8359], 'filled');
errorbar(us_pooled.reward.avg, us_pooled.punish.avg, -us_pooled.punish.SEM, us_pooled.punish.SEM,...
    -us_pooled.reward.SEM, us_pooled.reward.SEM, '.', 'Color', [0 0 0]);

% xlabel('Reward (\fontsize{20}\sigma\fontsize{16}-baseline)'); ylabel('Punish (\fontsize{20}\sigma\fontsize{16}-baseline)'); 
% xlabel('Reward (\sigma-baseline)'); ylabel('Punish (\sigma-baseline)'); 
ylim = get(gca, 'YLim'); xlim = get(gca, 'XLim');
set(gca, 'XTick', [0 1], 'YTick', [0 2], 'YLim', [-0.5 ylim(2)], 'XLim', [-0.5 xlim(2)]);
addOrginLines(gca, [0.7 0.7 0.7]);
xlabel('Reward'); ylabel('Punish');
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
    saveas(gcf, fullfile(figPath, [saveName '.fig']));
    saveas(gcf, fullfile(figPath, [saveName '.jpg']));   
end  

%% cs scatter plot
    % cs scatter plot
saveName = 'CuedOutcome_Cs_scatterPlot';
ensureFigure(saveName, 1); 
figSize = [1 1];

scatter(cs_pooled.low.avg,cs_pooled.high.avg, 36, [0 0 0], 'filled');
errorbar(cs_pooled.low.avg,cs_pooled.high.avg, -cs_pooled.high.SEM, cs_pooled.high.SEM,...
    -cs_pooled.low.SEM, cs_pooled.low.SEM, '.', 'Color', [0 0 0]);
% xlabel('Low Value (\fontsize{20}\sigma\fontsize{16}-baseline)'); ylabel('High Value (\fontsize{20}\sigma\fontsize{16}-baseline)');

set(gca, 'YLim', [0 1]); set(gca, 'XLim', [0 0.5]);
% setXYsameLimit;
set(gca, 'XTick', [0 0.5 1], 'YTick', [0 1]);
addUnityLine(gca, [0.7 0.7 0.7]);
xlabel('Low value'); ylabel('High value');
formatFigurePublish('size', figSize);
if saveOn
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
    saveas(gcf, fullfile(figPath, [saveName '.fig']));
    saveas(gcf, fullfile(figPath, [saveName '.jpg']));   
end  
    
    