% FrankenLNL_RewardPunish_poolAnimals

%{
goals:
1) Grand Averages per condition
2) Point estimates of prediction error, cue, outcome, by relevant condition
3) normalized valence coding
%}

DB = dbLoadExperiment('FrankenLNL_RewardPunish');

minRewardLickRate = 2; % at least n Hz licking (during us window)


%% Goal 1: Grand Averages
% data and associated descriptors for each data type, trial set combination
s3 = struct(...
    'data', [],...
    'animal', [],...
    'ch', []...
    );

% initialize different trial sets
s2 = struct(...
    'cuedReward', s3,...
    'uncuedReward', s3,...
    'omitReward', s3,...
    'cuedPuff', s3,...
    'uncuedPuff', s3,...
    'omitPuff', s3,...
    'cuedShock', s3,...
    'uncuedShock', s3,...
    'omitShock', s3...
    );

% s1
gAvg = struct(...
    'lickCue', s2,...
    'lick', s2,...
    'lickUs', s2,...    
    'phCue', s2,...
    'ph', s2,...
    'phUs', s2...
    );
    
for counter = 1:length(DB.animals)
    animal = DB.animals{counter}
%     if strcmp(animal, 'ACh_3')
%         continue;
%     end
    dbLoadAnimal(DB, animal);

    %% for each animal, assemble trial sets, do this manually to allow tweaking of trial makeup (exclude conditions such as late licking, no licking, select certain block numbers, etc.)
    %first column matches fields in gAvg, second column contains trial vector
    trialSets = {};
    
    % cuedReward
    trialSets{end+1, 1} = 'cuedReward';
    trialSets{end,2} = rewardTrials & Odor2Valve1Trials & (TE.licks_us.rate > minRewardLickRate);
    trialSets{end,3} = 3;
    % uncuedReward
    trialSets{end+1, 1} = 'uncuedReward';
    trialSets{end,2} = uncuedReward & (TE.licks_us.rate > minRewardLickRate);    
    trialSets{end,3} = 3;    
    % omitReward
    trialSets{end+1, 1} = 'omitReward';
    trialSets{end,2} = neutralTrials & Odor2Valve1Trials & (TE.licks_us.rate > minRewardLickRate);        
    trialSets{end,3} = 3;    
    % cuedPuff
    trialSets{end+1, 1} = 'cuedPuff';
    trialSets{end,2} = punishTrials & Odor2Valve2Trials;
    trialSets{end,3} = 3;    
    % uncuedPuff
    trialSets{end+1, 1} = 'uncuedPuff';
    trialSets{end,2} = uncuedPunish;
    trialSets{end,3} = 3;    
    % omitPuff
    trialSets{end+1, 1} = 'omitPuff';
    trialSets{end,2} = neutralTrials & Odor2Valve2Trials & (TE.BlockNumber == 3);
    trialSets{end,3} = 3;    
    % cuedShock
    trialSets{end+1, 1} = 'cuedShock';
    trialSets{end,2} = shockTrials & Odor2Valve2Trials;
    trialSets{end,3} = 4;    
    % uncuedShock
    trialSets{end+1, 1} = 'uncuedShock';
    trialSets{end,2} = uncuedShock;
    trialSets{end,3} = 4;     
    % omitShock
    trialSets{end+1, 1} = 'omitShock';
    trialSets{end,2} = neutralTrials & Odor2Valve2Trials & (TE.BlockNumber == 4);
    trialSets{end,3} = 4;     
    
    for tcounter = 1:size(trialSets, 1)
        trials = trialSets{tcounter, 2};
        label = trialSets{tcounter, 1};
        fullBlock = trialSets{tcounter, 3};
        % windows
        cueWindow = [zeros(size(TE.filename)) - 4 cellfun(@(x,y) y(end) - x(1), TE.Cue2, TE.Trace2)]; % allow for variable trace 2 duration, relative to cue onset
        window = [-7 4]; % relative to outcome
        usWindow = [0 4]; % relative to outcome

        % licking, cue        
        avgData = eventAverageFromTE(TE, trials, 'Port1In', 'zeroTimes', TE.Cue2, 'window', cueWindow);
        gAvg.lickCue.(label).data = expandVertCat(gAvg.lickCue.(label).data, avgData.Avg, 'left');
        gAvg.lickCue.(label).animal{end+1,1} = animal;
        gAvg.lickCue.(label).ch(end+1,1) = NaN;
        % licking, full (only relevant for block number = fullBlock)
        avgData = eventAverageFromTE(TE, trials & (TE.BlockNumber == fullBlock), 'Port1In', 'zeroTimes', TE.Us, 'window', window);
        gAvg.lick.(label).data = expandVertCat(gAvg.lick.(label).data, avgData.Avg, 'left');
        gAvg.lick.(label).animal{end+1,1} = animal;
        gAvg.lick.(label).ch(end+1,1) = NaN;                
        % licking, us
        avgData = eventAverageFromTE(TE, trials, 'Port1In', 'zeroTimes', TE.Us, 'window', usWindow);
        gAvg.lickUs.(label).data = expandVertCat(gAvg.lickUs.(label).data, avgData.Avg, 'left');
        gAvg.lickUs.(label).animal{end+1,1} = animal;
        gAvg.lickUs.(label).ch(end+1,1) = NaN;                
        
        
        for ch = 1:length(TE.Photometry.data)
            % photometry, cue        
            avgData = phAverageFromTE(TE, trials, ch, 'zeroTimes', TE.Cue2, 'window', cueWindow, 'FluorDataField', 'ZS');
            gAvg.phCue.(label).data = expandVertCat(gAvg.phCue.(label).data, avgData.Avg, 'left');
            gAvg.phCue.(label).animal{end+1,1} = animal;
            gAvg.phCue.(label).ch(end+1,1) = ch;
            % photometry, full (only relevant for block number = fullBlock)
            avgData = phAverageFromTE(TE, trials & (TE.BlockNumber == fullBlock), ch, 'zeroTimes', TE.Us, 'window', window, 'FluorDataField', 'ZS');
            gAvg.ph.(label).data = expandVertCat(gAvg.ph.(label).data, avgData.Avg, 'left');
            gAvg.ph.(label).animal{end+1,1} = animal;
            gAvg.ph.(label).ch(end+1,1) = ch;                
            % photometry, us
            avgData = phAverageFromTE(TE, trials, ch, 'zeroTimes', TE.Us, 'window', usWindow, 'FluorDataField', 'ZS');
            gAvg.phUs.(label).data = expandVertCat(gAvg.phUs.(label).data, avgData.Avg, 'left');
            gAvg.phUs.(label).animal{end+1,1} = animal;
            gAvg.phUs.(label).ch(end+1,1) = ch;      
        end
    end
end

%% normalize grand Average data to uncued reward us amplitude
gAvgNorm = gAvg;
% photometry reward Us size per fiber
normVectorPh = max(gAvg.phUs.uncuedReward.data(:,1:bpX2pnt(1, 20)), [], 2);
% lick reward Us size per fiber
% normVectorLick = nanmean(gAvg.lickUs.uncuedReward.data(:,1:bpX2pnt(1, 20)), 2);
for tcounter = 1:size(trialSets, 1)
    label = trialSets{tcounter, 1};
    gAvgNorm.phCue.(label).data = gAvg.phCue.(label).data ./ normVectorPh;
    gAvgNorm.ph.(label).data = gAvg.ph.(label).data ./ normVectorPh;
    gAvgNorm.phUs.(label).data = gAvg.phUs.(label).data ./ normVectorPh;
%     gAvgNorm.lickCue(label).data = gAvgNorm.lickCue(label).data ./ normVectorLick;
%     gAvgNorm.lick(label).data = gAvgNorm.lick(label).data ./ normVectorLick;
%     gAvgNorm.lickUs(label).data = gAvgNorm.lickUs(label).data ./ normVectorLick;
end
%% plot grand Averages
ensureFigure('raw', 1); plot(nanmean([gAvg.phCue.cuedReward.data gAvg.phUs.cuedReward.data])); hold on; plot(nanmean([gAvg.phCue.uncuedReward.data gAvg.phUs.uncuedReward.data]))
ensureFigure('norm', 1); plot(nanmean([gAvgNorm.phCue.cuedReward.data gAvgNorm.phUs.cuedReward.data])); hold on; plot(nanmean([gAvgNorm.phCue.uncuedReward.data gAvgNorm.phUs.uncuedReward.data]))


    