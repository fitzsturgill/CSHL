% script to pool data across experiments, collecting uncued Reward and Air
% Puff (Punishment) responses.

usWindow = [0 1]; % calculate mean of signal from 0 to n seconds post-Us
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
    RewardPunish_data(acounter + offset).Reward = response.data(rewardTrials & missTrials);
    RewardPunish_data(acounter + offset).Punish = response.data(punishTrials & CRTrials);
    RewardPunish_data(acounter + offset).animal = animal;
    RewardPunish_data(acounter + offset).experiment = 'LNL_reversals_punish';
end


%% save the data
    save(fullfile(savepath, 'RewardPunish_data.mat'), 'RewardPunish_data');
    disp(['*** Saved: ' fullfile(savepath, 'RewardPunish_data.mat')]);