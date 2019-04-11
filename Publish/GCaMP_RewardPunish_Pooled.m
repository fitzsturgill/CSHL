% script to pool data across experiments, collecting uncued Reward and Air
% Puff (Punishment) responses.

usWindow = [.1 0.6]; % calculate mean of signal from 0 to n seconds post-Us
savepath = 'Z:\SummaryAnalyses\ChAT_PE_Manuscript\GCaMP_RewardPunish_Pooled\';

CuedOutcome_animals = {'ChAT_26', 'ChAT_32', 'ChAT_34', 'ChAT_35', 'ChAT_37', 'ChAT_39', 'ChAT_42'};
SO_animals = {'ChAT_20', 'ChAT_21', 'ChAT_22'};
LNL_animals = {'DC_17', 'DC_20', 'DC_35', 'DC_36', 'DC_37', 'DC_40'};
totalAnimals = sum([length(CuedOutcome_animals) length(SO_animals) length(LNL_animals)]);
RewardPunish_data = struct(...
    'Reward', [],...
    'Punish', [],...
    'name', [],...
    'experiment', []...
    );

RewardPunish_data = repmat(RewardPunish_data, totalAnimals, 1);


% first CuedOutcome animals
DB = dbLoadExperiment('cuedOutcome');
for acounter = 1:length(CuedOutcome_animals)
    animal = CuedOutcome_animals{acounter};
    success = dbLoadAnimal(DB, animal); % load TE and trial lookups    
    response = bpCalcPeak_dFF(TE.Photometry, 1, usWindow, TE.Us, 'method', 'mean', 'phField', 'ZS');
    RewardPunish_data(acounter).Reward = response.data(filterTE(TE, 'trialType', 7, 'reject', 0));
    RewardPunish_data(acounter).Punish = response.data(filterTE(TE, 'trialType', 8, 'reject', 0));
    RewardPunish_data(acounter).animal = animal;
    RewardPunish_data(acounter).experiment = 'CuedOutcome_Odor_Complete';
end

%%
% next SO_RewardPunish animals
DB = dbLoadExperiment('SO_RewardPunish_odor');
offset = length(CuedOutcome_animals);
for acounter = 1:length(SO_animals)
    animal = SO_animals{acounter};
    success = dbLoadAnimal(DB, animal); % load TE and trial lookups    
    response = bpCalcPeak_dFF(TE.Photometry, 1, usWindow, TE.PostUsRecording, 'method', 'mean', 'phField', 'ZS');
    RewardPunish_data(acounter + offset).Reward = response.data(filterTE(TE, 'trialType', 5, 'reject', 0));
    RewardPunish_data(acounter + offset).Punish = response.data(filterTE(TE, 'trialType', 6, 'reject', 0));
    RewardPunish_data(acounter + offset).animal = animal;
    RewardPunish_data(acounter + offset).experiment = 'SO_RewardPunish_odor';
end

%% 
% last lickNoLick_odor_v2 with pavlovian_reversals_blocks (include air
% puff)
% I haven't created a database for this experiment yet....
basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\';
offset = sum([length(CuedOutcome_animals) length(SO_animals)]);
% clear counter; % get rid of variable in base workspace (created via experiment-specific "conditions" script, e.g. cuedOutcome_Conditions.m) so I can use "counter" as a name and don't have to do "acounter"

for acounter = 1:length(LNL_animals)
    animal = LNL_animals{acounter};
    load(fullfile(basepath, animal, 'TE.mat'));
    response = bpCalcPeak_dFF(TE.Photometry, 1, usWindow, TE.Us, 'method', 'mean', 'phField', 'ZS');
    LNL_conditions;
    if sum(uncuedReward)
        RewardPunish_data(acounter + offset).Reward = response.data(uncuedReward | (rewardTrials & missTrials)); %rewardTrials & missTrials
        RewardPunish_data(acounter + offset).Punish = response.data(uncuedPunish | (punishTrials & CRTrials)); %punishTrials & CRTrials
    else
        RewardPunish_data(acounter + offset).Reward = response.data(rewardTrials & missTrials); %rewardTrials & missTrials
        RewardPunish_data(acounter + offset).Punish = response.data(punishTrials & CRTrials); %punishTrials & CRTrials
    end
    RewardPunish_data(acounter + offset).animal = animal;
    RewardPunish_data(acounter + offset).experiment = 'LNL_reversals_punish';
end


%% save the data
    save(fullfile(savepath, 'RewardPunish_data.mat'), 'RewardPunish_data');
    disp(['*** Saved: ' fullfile(savepath, 'RewardPunish_data.mat')]);
    
    
%% bar graphs of reward vs punishment responses per animal with error bars
saveName = 'RewardPunish_barGraph_pooled';
ensureFigure(saveName, 1);

% compute mean and SEM for each animal, reward and punishment
clear counter; clear pooled;
nSubjects = length(RewardPunish_data);
pooled = struct(...
    'reward_mean', zeros(nSubjects, 1),...
    'reward_sem', zeros(nSubjects, 1),...
    'punish_mean', zeros(nSubjects, 1),...
    'punish_sem', zeros(nSubjects, 1),...
    'animal', [],...
    'experiment', []...    
    );
for counter = 1:nSubjects
    pooled.reward_mean(counter) = nanmean(RewardPunish_data(counter).Reward);
    pooled.reward_sem(counter) = nanstd(RewardPunish_data(counter).Reward) ./ sqrt(sum(isfinite(RewardPunish_data(counter).Reward)));
    sum(isfinite(RewardPunish_data(counter).Reward))
    sum(isfinite(RewardPunish_data(counter).Punish))
    pooled.punish_mean(counter) = nanmean(RewardPunish_data(counter).Punish);
    pooled.punish_sem(counter) = nanstd(RewardPunish_data(counter).Punish) ./ sqrt(sum(isfinite(RewardPunish_data(counter).Punish)));
    pooled.animal{counter} = RewardPunish_data(counter).animal;
    pooled.experiment{counter} = RewardPunish_data(counter).experiment;
end

axes;
% errorbar(pooled.reward_mean, pooled.punish_mean, pooled.punish_mean - pooled.punish_sem, pooled.punish_mean + pooled.punish_sem,...
%     pooled.reward_mean - pooled.reward_sem, pooled.reward_mean + pooled.reward_sem, 'or');
% scatter(pooled.reward_mean, pooled.punish_mean);
% set(gca, 'YScale', 'log', 'XScale', 'log');
bar((1:2:31) - 0.3, pooled.reward_mean, 'FaceColor', 'b', 'BarWidth', 0.25); hold on;
errorbar((1:2:31) - 0.3, pooled.reward_mean, pooled.reward_sem, '.k');
bar((1:2:31) + 0.3, pooled.punish_mean, 'FaceColor', 'r', 'BarWidth', 0.25); hold on;
errorbar((1:2:31) + 0.3, pooled.punish_mean, pooled.reward_sem, '.k');
% setXYsameLimit;
xlabel('Animal'); ylabel('Us (ZS)');
set(gca, 'XTickLabel', '');
% addUnityLine;

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
end    

saveName = 'RP_histograms_tiled';
fig = ensureFigure(saveName, 1);
params = struct();
params.cellmargin = [.05 .05 0.05 0.05];    
params.figmargin = [0.05 0.05 0.05 0.05];

hax = axesmatrix(4,4,1:nSubjects, params, fig);   % nSubjects = 16 here..
for counter = 1:nSubjects 
    axes(hax(counter)); hold on;
    histogram(RewardPunish_data(counter).Reward, -5:0.2:8, 'FaceColor', 'b', 'EdgeColor', 'none');
    histogram(RewardPunish_data(counter).Punish, -5:0.2:8, 'FaceColor', 'r', 'EdgeColor', 'none');
    textBox(RewardPunish_data(counter).animal);
end
% sameYScale(hax);

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
end    

% 4,7,10
saveName = 'RP_histograms_tiled_rewardOnly_examples';
fig = ensureFigure(saveName, 1);
params = struct();
params.cellmargin = [.05 .05 0.2 0.2];    
params.figmargin = [0.05 0.05 0.05 0.05];

hax = axesmatrix(1,3,1:3, params, fig);   % nSubjects = 16 here..
rewardOnly_examples = [4 7 10];
for counter = 1:length(rewardOnly_examples)
    axes(hax(counter)); hold on;
    ix = rewardOnly_examples(counter);
    histogram(RewardPunish_data(ix).Reward, -5:0.2:8, 'FaceColor', 'b', 'EdgeColor', 'none');
    histogram(RewardPunish_data(ix).Punish, -5:0.2:8, 'FaceColor', 'r', 'EdgeColor', 'none');
    textBox(RewardPunish_data(counter).animal, gca, [], 8);
end
axes(hax(2));
xlabel('Us response (mean ZS)');
% sameYScale(hax);
formatFigurePublish('size', [4 1.5]);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    export_fig(fullfile(savepath, saveName), '-eps');
end  




%% scatter plot of reward and punishment responses

saveName = 'RewardPunish_scatter_pooled';
ensureFigure(saveName, 1);

% compute mean and SEM for each animal, reward and punishment

axes;
% scatter(pooled.reward_mean, pooled.punish_mean, 10, mycolors('chat'), '.');
errorbar(pooled.reward_mean, pooled.punish_mean, -pooled.punish_sem, pooled.punish_sem,...
    -pooled.reward_sem, pooled.reward_sem, '.', 'Color', mycolors('chat'));
setXYsameLimit;
% ylims = get(gca, 'YLim');
% xlims = get(gca, 'XLim');
% set(gca, 'XLim', [-xlims(2) xlims(2)], 'YLim', [-ylims(2) ylims(2)]);
xlabel('Reward (ZS)'); ylabel('Punish (ZS)');
addOrginLines;
formatFigurePublish('size', [2 2]);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    export_fig(fullfile(savepath, saveName), '-eps');    
end    


