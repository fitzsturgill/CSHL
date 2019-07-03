%{
FrankenLNL_RewardPunish_poolAnimals
goals:
1) Grand Averages per condition
2) trial estimates of reward response
%}


DB = dbLoadExperiment('FrankenLNL_varyRewardSize');

saveOn = 1;
rewardLickBlankTime = 0.2;
minRewardLickCount = 2; % at least n licks in receipt of reward for reward conditions
lickUsWindow = [rewardLickBlankTime 1];
usWindow = [0.2 1];
grandAvgWindow = [-4 4];
baselineWindow = [-4 0]; % prior to Us
savepath = fullfile(DB.path, ['pooled' filesep]);
ensureDirectory(savepath);
% figsavepath = fullfile(DB.path, ['pooled' filesep 'figure']);
% ensureDirectory(figsavepath);
na = length(DB.animals);

% setup for grand averages
% data and associated descriptors for each data type, trial set combination
s3 = struct(...
    'data', [],...
    'animal', [],...
    'ch', []...
    );

% initialize different trial sets
s2 = struct(...
    'large', s3,...
    'medium', s3,...
    'small', s3,...
    'neutral', s3...
    );

% s1
gAvg = struct(...
    'lick', s2,...
    'ph', s2...
    );

% setup for us point estimates
s2 = struct(...
    'phPeakMean', zeros(na,2),...
    'phPeakMean_sem', zeros(na,2),...
    'phPeakPercentile', zeros(na, 2),...
    'phPeakPercentile_sem', zeros(na, 2),...
    'lickRate', zeros(na, 1),...
    'lickRate_sem', zeros(na, 1),...
    'n', zeros(na, 1)...
    );

us_pooled = struct(...
    'large', s2,...
    'medium', s2,...
    'small', s2,...
    'neutral', s2,...
    'animals', []...  
    );

for counter = 1:length(DB.animals)
    animal = DB.animals{counter};

    dbLoadAnimal(DB, animal);
    channels = TE.Photometry.settings.channels;
    TE.licks_us = countEventFromTE(TE, 'Port1In', lickUsWindow, TE.Us);
    TE.lick_baseline = countEventFromTE(TE, 'Port1In', baselineWindow, TE.Us);
    for channel = channels
        TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, TE.Us, 'method', 'percentile', 'percentile', 0.9, 'phField', 'ZS');
        TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, baselineWindow, TE.Us, 'method', 'mean', 'phField', 'ZS');
    end
    

    %% for each animal, assemble trial sets, do this manually to allow tweaking of trial makeup, min lick counts, etc.
    %first column matches fields in gAvg, second column contains trial vector
    trialSets = {};
    
    % large reward
    trialSets{end+1, 1} = 'large';
    trialSets{end,2} = largeRewardTrials & (TE.licks_us.count > minRewardLickCount) & (TE.BlockNumber == 10);
    
    % medium reward
    trialSets{end+1, 1} = 'medium';
    trialSets{end,2} = mediumRewardTrials & (TE.licks_us.count > minRewardLickCount) & (TE.BlockNumber == 10);
    
    % large reward
    trialSets{end+1, 1} = 'small';
    trialSets{end,2} = smallRewardTrials & (TE.licks_us.count > minRewardLickCount) & (TE.BlockNumber == 10);    
    
    % neutral
    trialSets{end+1, 1} = 'neutral';
    trialSets{end,2} = neutralTrials & (TE.BlockNumber == 10);        
    
    for tcounter = 1:size(trialSets, 1)
        trials = trialSets{tcounter, 2};        
        label = trialSets{tcounter, 1};
        
        % licking, us
        avgData = eventAverageFromTE(TE, trials, 'Port1In', 'zeroTimes', TE.Us, 'window', grandAvgWindow);
        gAvg.lick.(label).data = expandVertCat(gAvg.lick.(label).data, avgData.Avg, 'left');
        gAvg.lick.(label).animal{end+1,1} = animal;
        gAvg.lick.(label).ch(end+1,1) = NaN;
        
%         % setup for us point estimates

        nt = sum(trials);
        us_pooled.(label).n(counter) = nt;
        % us_pooled
        us_pooled.(label).lickRate(counter) = nanmean(TE.licks_us.rate(trials));
        us_pooled.(label).lickRate_sem(counter) = nanstd(TE.licks_us.rate(trials)) / sqrt(nt);
        
       for ch = 1:length(TE.Photometry.data)
            % photometry, full (only relevant for block number = fullBlock)
            avgData = phAverageFromTE(TE, trials, ch, 'zeroTimes', TE.Us, 'window', grandAvgWindow, 'FluorDataField', 'ZS');
            gAvg.ph.(label).data = expandVertCat(gAvg.ph.(label).data, avgData.Avg, 'left');
            gAvg.ph.(label).animal{end+1,1} = animal;
            gAvg.ph.(label).ch(end+1,1) = ch;                

            % us_pooled
            us_pooled.(label).phPeakMean(counter, ch) = nanmean(TE.phPeakMean_us(ch).data(trials));
            us_pooled.(label).phPeakMean_sem(counter, ch) = nanstd(TE.phPeakMean_us(ch).data(trials)) / sqrt(nt);

            us_pooled.(label).phPeakPercentile(counter, ch) = nanmean(TE.phPeakPercentile_us(ch).data(trials));
            us_pooled.(label).phPeakPercentile_sem(counter, ch) = nanstd(TE.phPeakPercentile_us(ch).data(trials)) / sqrt(nt);            
       end
    
    end       
end

%% compute averages, SEM, 

for tcounter = 1:size(trialSets, 1)
    label = trialSets{tcounter, 1}; 
    gAvg.ph.(label).Avg = nanmean(gAvg.ph.(label).data);
    gAvg.ph.(label).SEM = nanstd(gAvg.ph.(label).data, 0, 1) ./ sum(isfinite(gAvg.ph.(label).data), 1);
    gAvg.ph.(label).xData = (0:(size(gAvg.ph.(label).data, 2) - 1)) * 1/20 - (size(gAvg.ph.(label).data, 2)/20 - 4);    
end

save(fullfile(savepath, 'grandAverages.mat'), 'gAvg');
disp(['*** saving: ' fullfile(savepath, 'grandAverages.mat') ' ***']);

save(fullfile(savepath, 'us_pooled.mat'), 'us_pooled');
disp(['*** saving: ' fullfile(savepath, 'us_pooled.mat') ' ***']);

