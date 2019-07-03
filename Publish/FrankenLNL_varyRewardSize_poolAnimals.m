%{
FrankenLNL_RewardPunish_poolAnimals
goals:
1) Grand Averages per condition
2) trial estimates of reward response
%}


DB = dbLoadExperiment('FrankenLNL_varyRewardSize');

saveOn = 1;
minRewardLickLatency = 2; % at least n Hz licking (during us window for reward trials)
rewardLickBlankTime = 0.2;
savepath = fullfile(DB.path, ['pooled' filesep]);
ensureDirectory(savepath);
figsavepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(figsavepath);

%% Goal 1: Grand Averages
% data and associated descriptors for each data type, trial set combination
s3 = struct(...
    'data', [],...
    'animal', [],...
    'ch', []...
    );

% initialize different trial sets
s2 = struct(...
    'largeReward', s3,...
    'mediumReward', s3,...
    'smallReward', s3,...
    'neutral', s3...
    );

% s1
gAvg = struct(...
    'lickBaseline', s2,...
    'lickUs', s2,...
    'phBaseline', s2,...
    'phUs', s2...
    );


for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
%     if strcmp(animal, 'ACh_3')
%         continue;
%     end
    dbLoadAnimal(DB, animal);

    %% for each animal, assemble trial sets, do this manually to allow tweaking of trial makeup (exclude conditions such as late licking, no licking, select certain block numbers, etc.)
    %first column matches fields in gAvg, second column contains trial vector
    trialSets = {};
    
    % cuedReward
    trialSets{end+1, 1} = 'largeReward';
    trialSets{end,2} = rewardTrials & Odor2Valve1Trials & (TE.licks_us.rate > minRewardLickRate) & (TE.licks_cs.rate > minCueLickRate) & ismember(TE.BlockNumber, [2 3]); % exclude shock or extinction days
    trialSets{end,3} = 3;
    
    
end
