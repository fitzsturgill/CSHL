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

%% compute averages, SEM, also normalize grand Average data to uncued reward us amplitude
gAvgNorm = gAvg;
% photometry reward Us size per fiber
normVectorPh = max(gAvg.phUs.uncuedReward.data(:,1:bpX2pnt(1, 20)), [], 2);
% lick reward Us size per fiber
% normVectorLick = nanmean(gAvg.lickUs.uncuedReward.data(:,1:bpX2pnt(1, 20)), 2);
for tcounter = 1:size(trialSets, 1)
    label = trialSets{tcounter, 1};
    gAvg.phCue.(label).Avg = nanmean(gAvg.phCue.(label).data);
    gAvg.phCue.(label).SEM = nanstd(gAvg.phCue.(label).data, 0, 1) ./ sum(isfinite(gAvg.phCue.(label).data), 1);
    gAvg.phCue.(label).xData = (0:(size(gAvg.phCue.(label).data, 2) - 1)) * 1/20 - 7;
    gAvg.ph.(label).Avg = nanmean(gAvg.ph.(label).data);
    gAvg.ph.(label).SEM = nanstd(gAvg.ph.(label).data, 0, 1) ./ sum(isfinite(gAvg.ph.(label).data), 1);
    gAvg.ph.(label).xData = (0:(size(gAvg.phCue.(label).data, 2) - 1)) * 1/20 - (size(gAvg.phCue.(label).data, 2)/20 - 4);
    gAvg.phUs.(label).Avg = nanmean(gAvg.phUs.(label).data);
    gAvg.phUs.(label).SEM = nanstd(gAvg.phUs.(label).data, 0, 1) ./ sum(isfinite(gAvg.phUs.(label).data), 1);    
    gAvg.phUs.(label).xData = (0:(size(gAvg.phUs.(label).data, 2) - 1)) * 1/20;
    
    gAvgNorm.phCue.(label).data = gAvg.phCue.(label).data ./ normVectorPh;
    gAvgNorm.ph.(label).data = gAvg.ph.(label).data ./ normVectorPh;
    gAvgNorm.phUs.(label).data = gAvg.phUs.(label).data ./ normVectorPh;
    
    gAvgNorm.phCue.(label).Avg = nanmean(gAvgNorm.phCue.(label).data);
    gAvgNorm.phCue.(label).SEM = nanstd(gAvgNorm.phCue.(label).data, 0, 1) ./ sum(isfinite(gAvgNorm.phCue.(label).data), 1);
    gAvgNorm.phCue.(label).xData = (0:(size(gAvgNorm.phCue.(label).data, 2) - 1)) * 1/20 - 7;
    gAvgNorm.ph.(label).Avg = nanmean(gAvgNorm.ph.(label).data);
    gAvgNorm.ph.(label).SEM = nanstd(gAvgNorm.ph.(label).data, 0, 1) ./ sum(isfinite(gAvgNorm.ph.(label).data), 1);
    gAvgNorm.ph.(label).xData = (0:(size(gAvgNorm.phCue.(label).data, 2) - 1)) * 1/20 - (size(gAvgNorm.phCue.(label).data, 2)/20 - 4);
    gAvgNorm.phUs.(label).Avg = nanmean(gAvgNorm.phUs.(label).data);
    gAvgNorm.phUs.(label).SEM = nanstd(gAvgNorm.phUs.(label).data, 0, 1) ./ sum(isfinite(gAvgNorm.phUs.(label).data), 1);    
    gAvgNorm.phUs.(label).xData = (0:(size(gAvgNorm.phUs.(label).data, 2) - 1)) * 1/20;
    
%     gAvgNorm.lickCue(label).data = gAvgNorm.lickCue(label).data ./ normVectorLick;
%     gAvgNorm.lick(label).data = gAvgNorm.lick(label).data ./ normVectorLick;
%     gAvgNorm.lickUs(label).data = gAvgNorm.lickUs(label).data ./ normVectorLick;
end
return;
%% plot grand Averages
ensureFigure('grandAverages', 1); 
% xData = [0:(size(gAvg.phCue.cuedReward.data, 2) - 1) * 1/20 - 4 ...
%     0:(size(gAvg.phUs.cuedReward.data, 2) - 1) * 1/20];
subplot(2,2,1); %plot(nanmean([gAvg.phCue.cuedReward.data gAvg.phUs.cuedReward.data])); hold on; plot(nanmean([gAvg.phCue.uncuedReward.data gAvg.phUs.uncuedReward.data])); 
boundedline([gAvg.phCue.cuedReward.xData gAvg.phUs.cuedReward.xData],...
    [[gAvg.phCue.cuedReward.Avg'; gAvg.phUs.cuedReward.Avg'] [gAvg.phCue.uncuedReward.Avg'; gAvg.phUs.uncuedReward.Avg']],...
    permute([[gAvg.phCue.cuedReward.SEM gAvg.phUs.cuedReward.SEM] ; [gAvg.phCue.cuedReward.SEM gAvg.phUs.cuedReward.SEM]], [2 3 1]));
title('appetitive');
subplot(2,2,2); plot(nanmean([gAvgNorm.phCue.cuedReward.data gAvgNorm.phUs.cuedReward.data])); hold on; plot(nanmean([gAvgNorm.phCue.uncuedReward.data gAvgNorm.phUs.uncuedReward.data]))
title('appetitive norm.');
subplot(2,2,3); plot(nanmean([gAvgNorm.phCue.cuedPuff.data gAvgNorm.phUs.cuedPuff.data])); hold on; plot(nanmean([gAvgNorm.phCue.uncuedPuff.data gAvgNorm.phUs.uncuedPuff.data]))
title('normPunish');
subplot(2,2,4); plot(nanmean([gAvgNorm.phCue.cuedShock.data gAvgNorm.phUs.cuedShock.data])); hold on; plot(nanmean([gAvgNorm.phCue.uncuedShock.data gAvgNorm.phUs.uncuedShock.data]))
title('normShock');


%% 
boundedline(xData(:), [nanmean(pupData(trialsByType{3}, :)); nanmean(pupData(trialsByType{4} & TE.BlockNumber == 4, :)); nanmean(pupData(trialsByType{6}, :))]',...
    permute([nanstd(pupData(trialsByType{3}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{3},:)), 1));...
    nanstd(pupData(trialsByType{4} & TE.BlockNumber == 4, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{4} & TE.BlockNumber == 4,:)), 1));...
    nanstd(pupData(trialsByType{6}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{6},:)), 1))], [2 3 1]),...
    'cmap', [1 0 0; 0 0 0; 1 0 1]);
    