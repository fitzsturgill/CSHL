% FrankenLNL_RewardPunish_poolAnimals

%{
goals:
1) Grand Averages per condition
2) Latency, reliability and jitter of Us responses
2) Point estimates of prediction error, cue, outcome, by relevant condition
3) normalized valence coding
%}

DB = dbLoadExperiment('FrankenLNL_RewardPunish');
saveOn = 1;
minRewardLickRate = 2; % at least n Hz licking (during us window for reward trials)
minCueLickRate = 1; % at least n Hz licking (during cs window for reward cued reward trials)
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
    'lickDelay', s2,...
    'lick', s2,...
    'lickUs', s2,...    
    'phCue', s2,...
    'phDelay', s2,...
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
    trialSets{end,2} = rewardTrials & Odor2Valve1Trials & (TE.licks_us.rate > minRewardLickRate) & (TE.licks_cs.rate > minCueLickRate) & ismember(TE.BlockNumber, [2 3]); % exclude shock or extinction days
    trialSets{end,3} = 3;
    % uncuedReward
    trialSets{end+1, 1} = 'uncuedReward';
    trialSets{end,2} = uncuedReward & (TE.licks_us.rate > minRewardLickRate) & ismember(TE.BlockNumber, [2 3]);    
    trialSets{end,3} = 3;    
    % omitReward
    trialSets{end+1, 1} = 'omitReward';
    trialSets{end,2} = neutralTrials & Odor2Valve1Trials & ismember(TE.BlockNumber, [2 3]);        
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
        % baseline, cue, and delay period: allow for variable trace 2 duration, relative to cue onset        
        cueWindow = [zeros(size(TE.filename)) - 4 cellfun(@(x,y) y(end) - x(1), TE.Cue2, TE.Trace2)]; 
        % just the delay period, relative to outcome
        delayWindow = [zeros(size(TE.filename)) - cellfun(@(x) x(end) - x(1), TE.Trace2) zeros(size(TE.filename))];
        window = [-7 4]; % relative to outcome
        usWindow = [0 4]; % relative to outcome

        % licking, cue        
        avgData = eventAverageFromTE(TE, trials, 'Port1In', 'zeroTimes', TE.Cue2, 'window', cueWindow);
        gAvg.lickCue.(label).data = expandVertCat(gAvg.lickCue.(label).data, avgData.Avg, 'left');
        gAvg.lickCue.(label).animal{end+1,1} = animal;
        gAvg.lickCue.(label).ch(end+1,1) = NaN;
        % licking, delay, aligned at right        
        avgData = eventAverageFromTE(TE, trials, 'Port1In', 'zeroTimes', TE.Us, 'window', delayWindow);
        gAvg.lickDelay.(label).data = expandVertCat(gAvg.lickDelay.(label).data, avgData.Avg, 'right');
        gAvg.lickDelay.(label).animal{end+1,1} = animal;
        gAvg.lickDelay.(label).ch(end+1,1) = NaN;        
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
            % photometry, delay        
            avgData = phAverageFromTE(TE, trials, ch, 'zeroTimes', TE.Us, 'window', delayWindow, 'FluorDataField', 'ZS');
            gAvg.phDelay.(label).data = expandVertCat(gAvg.phDelay.(label).data, avgData.Avg, 'right');
            gAvg.phDelay.(label).animal{end+1,1} = animal;
            gAvg.phDelay.(label).ch(end+1,1) = ch;            
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
    gAvg.phDelay.(label).Avg = nanmean(gAvg.phDelay.(label).data);
    gAvg.phDelay.(label).SEM = nanstd(gAvg.phDelay.(label).data, 0, 1) ./ sum(isfinite(gAvg.phDelay.(label).data), 1);
    gAvg.phDelay.(label).xData = (0:(size(gAvg.phDelay.(label).data, 2) - 1)) * 1/20 - 2;    
    gAvg.ph.(label).Avg = nanmean(gAvg.ph.(label).data);
    gAvg.ph.(label).SEM = nanstd(gAvg.ph.(label).data, 0, 1) ./ sum(isfinite(gAvg.ph.(label).data), 1);
    gAvg.ph.(label).xData = (0:(size(gAvg.phCue.(label).data, 2) - 1)) * 1/20 - (size(gAvg.phCue.(label).data, 2)/20 - 4);
    gAvg.phUs.(label).Avg = nanmean(gAvg.phUs.(label).data);
    gAvg.phUs.(label).SEM = nanstd(gAvg.phUs.(label).data, 0, 1) ./ sum(isfinite(gAvg.phUs.(label).data), 1);    
    gAvg.phUs.(label).xData = (0:(size(gAvg.phUs.(label).data, 2) - 1)) * 1/20;
    
    gAvgNorm.phCue.(label).data = gAvg.phCue.(label).data ./ normVectorPh;
    gAvgNorm.phDelay.(label).data = gAvg.phDelay.(label).data ./ normVectorPh;
    gAvgNorm.ph.(label).data = gAvg.ph.(label).data ./ normVectorPh;
    gAvgNorm.phUs.(label).data = gAvg.phUs.(label).data ./ normVectorPh;
    
    gAvgNorm.phCue.(label).Avg = nanmean(gAvgNorm.phCue.(label).data);
    gAvgNorm.phCue.(label).SEM = nanstd(gAvgNorm.phCue.(label).data, 0, 1) ./ sum(isfinite(gAvgNorm.phCue.(label).data), 1);
    gAvgNorm.phCue.(label).xData = (0:(size(gAvgNorm.phCue.(label).data, 2) - 1)) * 1/20 - 7;
    gAvgNorm.phDelay.(label).Avg = nanmean(gAvg.phDelay.(label).data);
    gAvgNorm.phDelay.(label).SEM = nanstd(gAvg.phDelay.(label).data, 0, 1) ./ sum(isfinite(gAvg.phDelay.(label).data), 1);
    gAvgNorm.phDelay.(label).xData = (0:(size(gAvg.phDelay.(label).data, 2) - 1)) * 1/20 - 2;        
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

save(fullfile(savepath, 'grandAverages.mat'), 'gAvg');
disp(['*** saving: ' fullfile(savepath, 'grandAverages.mat') ' ***']);
save(fullfile(savepath, 'grandAveragesNorm.mat'), 'gAvgNorm');
disp(['*** saving: ' fullfile(savepath, 'grandAveragesNorm.mat') ' ***']);


%%
savepath = fullfile(DB.path, ['pooled' filesep]);
load(fullfile(savepath, 'grandAverages.mat'), 'gAvg');
disp(['*** loading: ' fullfile(savepath, 'grandAverages.mat') ' ***']);
load(fullfile(savepath, 'grandAveragesNorm.mat'), 'gAvgNorm');
disp(['*** loading: ' fullfile(savepath, 'grandAveragesNorm.mat') ' ***']);


%% plot grand Averages


save(fullfile(savepath, 'grandAverages.mat'), 'gAvg');
disp(['*** saving: ' fullfile(savepath, 'grandAverages.mat') ' ***']);
save(fullfile(savepath, 'grandAveragesNorm.mat'), 'gAvgNorm');
disp(['*** saving: ' fullfile(savepath, 'grandAveragesNorm.mat') ' ***']);


saveName = 'grandAverage_appetitive_simple';
ensureFigure(saveName, 1); 
linecolors = [0 0 1; 0 0 0; 0 1 1];     
window = [-5 3];
axes; hold on;
xData = [gAvg.phCue.cuedReward.xData gAvg.phUs.cuedReward.xData;...
    gAvg.phCue.omitReward.xData gAvg.phUs.omitReward.xData;...
    gAvg.phCue.uncuedReward.xData gAvg.phUs.uncuedReward.xData]';
yData = [gAvg.phCue.cuedReward.Avg gAvg.phUs.cuedReward.Avg;...
    gAvg.phCue.omitReward.Avg gAvg.phUs.omitReward.Avg;
    gAvg.phCue.uncuedReward.Avg gAvg.phUs.uncuedReward.Avg]';
bData = permute([gAvg.phCue.cuedReward.SEM gAvg.phUs.cuedReward.SEM;...
    gAvg.phCue.omitReward.SEM gAvg.phUs.omitReward.SEM;...
    gAvg.phCue.uncuedReward.SEM gAvg.phUs.uncuedReward.SEM], [2 3 1]);
[hl, hp] = boundedline(xData, yData, bData, 'cmap', linecolors);
set(gca, 'XLim', window);
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% legend(hl, {'cued', 'omit', 'uncued'}, 'Location', 'best'); legend('boxoff');
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);
xlabel('Time from reinforcement (s)');

formatFigurePublish('size', [2 1]);

if saveOn 
    export_fig(fullfile(figsavepath, saveName), '-eps');
end


%% make avg rasters
% STUB

saveName = 'uncuedUs_avgRasters_rewNorm';
ensureFigure(saveName, 1);
clim = [-1 1];
xData = [-2 4]; yData = [1 length(DB.animals) + 0.5];
subplot(1,3,1); imagesc(xData, yData, [gAvgNorm.phDelay.uncuedReward.data gAvgNorm.phUs.uncuedReward.data], clim); title('reward'); ylabel('Mouse #');
subplot(1,3,2); imagesc(xData, yData, [gAvgNorm.phDelay.uncuedPuff.data gAvgNorm.phUs.uncuedPuff.data], clim); title('air puff'); set(gca, 'YTickLabel', {});
xlabel('Time from reinforcment (s)');
subplot(1,3,3); imagesc(xData, yData, [gAvgNorm.phDelay.uncuedShock.data gAvgNorm.phUs.uncuedShock.data], clim); title('shock'); set(gca, 'YTickLabel', {});
formatFigurePublish('size', [2 1]);
if saveOn 
    export_fig(fullfile(figsavepath, saveName), '-eps');
end

saveName = 'uncuedUs_avgRasters_rew';
ensureFigure(saveName, 1);
clim = [-1 5];
subplot(1,3,1); imagesc(xData, yData, [gAvg.phDelay.uncuedReward.data gAvg.phUs.uncuedReward.data], clim); title('reward'); ylabel('Mouse #');
subplot(1,3,2); imagesc(xData, yData, [gAvg.phDelay.uncuedPuff.data gAvg.phUs.uncuedPuff.data], clim); title('air puff'); set(gca, 'YTickLabel', {});
xlabel('Time from reinforcment (s)');
subplot(1,3,3); imagesc(xData, yData, [gAvg.phDelay.uncuedShock.data gAvg.phUs.uncuedShock.data], clim); title('shock'); set(gca, 'YTickLabel', {});
formatFigurePublish('size', [2 1]);
if saveOn 
    export_fig(fullfile(figsavepath, saveName), '-eps');
end

%% plot the averages for each animal

saveName = 'uncuedUs_avgs_rew';
ensureFigure(saveName, 1);
nAnimals = length(DB.animals);
nRows = ceil(sqrt(nAnimals));
nColumns = ceil(nAnimals/nRows);
xData = linspace(-2, 4, size([gAvg.phDelay.uncuedReward.data gAvg.phUs.uncuedReward.data], 2));
for counter = 1:nAnimals
    subplot(nRows, nColumns, counter);
    plot(xData, [gAvg.phDelay.uncuedReward.data(counter * 2 - 1:counter*2,:) gAvg.phUs.uncuedReward.data(counter * 2 - 1:counter*2,:)]');
end

saveName = 'uncuedUs_avgs_puff';
ensureFigure(saveName, 1);
for counter = 1:nAnimals
    subplot(nRows, nColumns, counter);
    plot(xData, [gAvg.phDelay.uncuedPuff.data(counter * 2 - 1:counter*2,:) gAvg.phUs.uncuedPuff.data(counter * 2 - 1:counter*2,:)]');
end

saveName = 'uncuedUs_avgs_shock';
ensureFigure(saveName, 1);
for counter = 1:nAnimals
    subplot(nRows, nColumns, counter);
    plot(xData, [gAvg.phDelay.uncuedShock.data(counter * 2 - 1:counter*2,:) gAvg.phUs.uncuedShock.data(counter * 2 - 1:counter*2,:)]');
end

%% calculate latency, jitter, and reliability of uncued Us responses

minRewardLickRate = 2;
bl_window = [-1 0];
response_window = [0 3]; % must start at 0
endTauOff = 0.8; % use up to n seconds after us delivery to fit tau off
PhotometryField = 'PhotometryHF';
na = length(DB.animals);

s2 = struct(...
    'latency', zeros(na, 2),...
    'latencyY', zeros(na, 2),... % 1/2 max value
    'avg_delta', zeros(na, 2),... % delta responsee from baseline
    'avgMax', NaN(na, 2),... % maximum (or minimum) value of response
    'tauOff', NaN(na, 2),...
    'tauObject', [],...
    'tauData', [],...
    'tauDataX', [],...
    'avgData', [],...
    'avgDataX', [],...
    'jitter_std', zeros(na, 2),...
    'jitter_avg', zeros(na, 2),...
    'jitter_tt50', [],... % cell array to hold distributions
    'jitter_delta', []... % cell array to hold distributions
    );

us_pooled = struct(...
    'rew', s2,...
    'puff', s2...
    );

trialSets = {'rew', 'puff', 'shock'};
rewAvgDelta = NaN(length(DB.animals), 2); % hold average rew delta ZS of uncued reward trials
puffAvgDelta = NaN(length(DB.animals), 2); % hold average rew delta ZS of uncued punish trials
for counter = 1:length(DB.animals)
    animal = DB.animals{counter}
    dbLoadAnimal(DB, animal);
    Fs = TE.(PhotometryField).sampleRate;
    for tscounter = 1:length(trialSets)
        trialSet = trialSets{tscounter};
        switch trialSet
            case 'rew'
                theseTrials = find(rewardTrials & uncuedTrials & (TE.licks_us.rate > minRewardLickRate));    
            case 'puff'
                theseTrials = find(punishTrials & uncuedTrials);   
            case 'shock'
                theseTrials = find(shockTrials & uncuedTrials);
        end
        [tt50, deltaMax] = deal(NaN(length(theseTrials), 2));
        for channel = 1:2
            [responseData, xData] = phAlignedWindow(TE, theseTrials, channel, 'zeroTimes', TE.Us, 'window', response_window, 'FluorDataField', 'ZS', 'PhotometryField', PhotometryField);
            [blData, blx] = phAlignedWindow(TE, theseTrials, channel, 'zeroTimes', TE.Us, 'window', bl_window, 'FluorDataField', 'ZS', 'PhotometryField', PhotometryField);
            bl = mean(blData, 2);
            response_avg = nanmean(responseData);
            % check if response is positive or negative
            if (max(response_avg) - mean(bl)) > (mean(bl) - min(response_avg))
                [peakVal, iAvg] = max(response_avg);
                avgDelta = peakVal - mean(bl);
                direction = 1;
                [M, I] = max(responseData, [], 2);                
                
                % fit exponential
                fitData = response_avg(iAvg:min(bpX2pnt(xData(iAvg) + endTauOff, Fs), length(response_avg)));
                fitXData = xData(iAvg:min(bpX2pnt(xData(iAvg) + endTauOff, Fs), length(response_avg)));            
                model = 'a + b * exp(-1/c*x)';
                fo = fitoptions('Method', 'NonlinearLeastSquares',...
                    'Upper', [min(response_avg), Inf, Inf],...
                    'Lower', [min(response_avg), 0, 0],...
                    'StartPoint', [min(response_avg), range(response_avg), 1]...
                    );
                ft = fittype(model, 'options', fo);
                [fitobject, gof, output] = fit(fitXData(:) - fitXData(1), fitData(:), ft, fo);
                us_pooled.(trialSet).tauOff(counter,channel) = fitobject.c;
                us_pooled.(trialSet).tauObject{counter,channel} = fitobject;
                us_pooled.(trialSet).tauData{counter,channel} = fitobject.a + fitobject.b * exp(-1/fitobject.c * (fitXData - fitXData(1)));
                us_pooled.(trialSet).tauDataX{counter,channel} = fitXData;                
            else
                [peakVal, iAvg] = min(response_avg);
                avgDelta = peakVal - mean(bl);
                direction = -1;
                if strcmp(trialSet, 'rew')
                    disp([animal 'wtf']);
                end
                [M, I] = min(responseData, [], 2);                
            end
            
            us_pooled.(trialSet).avgDataX{counter,channel} = xData;
            us_pooled.(trialSet).avgData{counter,channel} = response_avg;
                                                                           
            deltaMax(:,channel) = M - bl;
            level50 = bl + deltaMax(:,channel) * 0.5;
            
            [ind, to] = crossing(response_avg, xData, mean(bl) + avgDelta/2);
            us_pooled.(trialSet).latency(counter, channel) = to(1);
            us_pooled.(trialSet).latencyY(counter, channel) = avgDelta/2;
            us_pooled.(trialSet).avg_delta(counter,channel) = avgDelta; % can be negative or positive
            
            % find time to 50% per trial and magnitude of trial response
            thresh = responseData >= level50;
            crossings = [false(size(thresh,1), 1) thresh(:,2:end) - thresh(:,1:end-1)];
            crossings = crossings > 0; % only want upward crossings
            [row, col] = find(crossings);
            pointMatrix = NaN(size(crossings));
            pointMatrix(sub2ind(size(crossings), row,col)) = col;
            valid = isfinite(min(pointMatrix,[],2)); % only include trials where crossing is found
            theseCrossings = min(pointMatrix(valid,:),[],2);
            responseDataValid = responseData(valid,:);
            x2 = xData(theseCrossings)';
            x1 = xData(theseCrossings - 1)';
            y2 = responseDataValid(sub2ind(size(responseDataValid), (1:length(theseCrossings))', theseCrossings));
            y1 = responseDataValid(sub2ind(size(responseDataValid), (1:length(theseCrossings))', theseCrossings - 1)); 
            m = (y2 - y1) ./ (x2 - x1);
            x = (level50(valid) - y1 + m .* x1) ./ m;
            tt50(valid,channel) = x;       
    %         ensureFigure('test', 1);
    %         trials = randperm(length(theseTrials), 16);
    %         for counter = 1:16
    %             trial = trials(counter);
    %             subplot(4,4,counter); hold on;
    %             plot(xData, responseData(trial,:), '-*');
    %             plot(blx, blData(trial,:));
    %             scatter(tt50(trial, 2), level50(trial));
    %             plot(get(gca, 'XLim'), [level50(trial) level50(trial)], '--');
    %         end
        end       
        us_pooled.(trialSet).jitter_tt50{counter} = tt50;
        us_pooled.(trialSet).jitter_std(counter, :) = nanstd(tt50);
        us_pooled.(trialSet).jitter_avg(counter, :) = nanmean(tt50);
        us_pooled.(trialSet).jitter_delta{counter} = deltaMax;    
    end
end

%% plot scatter plots for reward and air puff of tt50 vs delta
minDelta = 2; % minimum deltaF in baseline standard deviation units (ZS)
saveName = 'jitter_tt50_vs_delta';
ensureFigure(saveName, 1);

rew_tt50_all = [];
puff_tt50_all = [];
shock_tt50_all = [];
for counter = 1:length(us_pooled.rew.jitter_tt50)
    for channel = 1:2
        if us_pooled.rew.avg_delta(counter, channel) > minDelta
            subplot(2,3,1); hold on;
            scatter(us_pooled.rew.jitter_tt50{counter}(:,channel), us_pooled.rew.jitter_delta{counter}(:,channel), '.');
            rew_tt50_all = [rew_tt50_all; us_pooled.rew.jitter_tt50{counter}(:,channel)];
        end
        if us_pooled.puff.avg_delta(counter, channel) > minDelta
            subplot(2,3,2); hold on;
            scatter(us_pooled.puff.jitter_tt50{counter}(:,channel), us_pooled.puff.jitter_delta{counter}(:,channel), '.');    
            puff_tt50_all = [puff_tt50_all; us_pooled.puff.jitter_tt50{counter}(:,channel)];
        end
        if us_pooled.shock.avg_delta(counter, channel) > minDelta
            subplot(2,3,3); hold on;
            scatter(us_pooled.shock.jitter_tt50{counter}(:,channel), us_pooled.shock.jitter_delta{counter}(:,channel), '.');    
            shock_tt50_all = [shock_tt50_all; us_pooled.shock.jitter_tt50{counter}(:,channel)];
        end        
    end
end
subplot(2,3,1); title('Reward'); ylabel('amplitude (ZS)');
subplot(2,3,2); title('Air Puff'); 
subplot(2,3,3); title('Shock'); 

binSize = 0.05;
bins = [0:binSize:2];
subplot(2,3,4); hold on; histogram(rew_tt50_all, bins); xlabel('latency (s)'); 
[N, edges] = histcounts(rew_tt50_all, bins); [~, ix] = max(N); rewMode = edges(ix) + binSize/2;
plot([rewMode rewMode], get(gca, 'YLim'), '--r'); 
set(gca, 'XLim', [0 max(get(gca, 'XLim'))]);
textBox(sprintf('mode=%.4g(s)', rewMode));

subplot(2,3,5); hold on; histogram(puff_tt50_all, bins); xlabel('latency (s)'); 
[N, edges] = histcounts(puff_tt50_all, bins); [~, ix] = max(N); puffMode = edges(ix) + binSize/2;
plot([puffMode puffMode], get(gca, 'YLim'), '--r'); 
set(gca, 'XLim', [0 max(get(gca, 'XLim'))]);
textBox(sprintf('mode=%.4g(s)', puffMode));

subplot(2,3,6); hold on; histogram(shock_tt50_all, bins); xlabel('latency (s)'); 
[N, edges] = histcounts(shock_tt50_all, bins); [~, ix] = max(N); shockMode = edges(ix) + binSize/2;
plot([shockMode shockMode], get(gca, 'YLim'), '--r'); 
set(gca, 'XLim', [0 max(get(gca, 'XLim'))]);
textBox(sprintf('mode=%.4g(s)', shockMode));

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
end        
%% plot averages with tauOff annotated




nAnimals = length(DB.animals);
nRows = ceil(sqrt(nAnimals));
nColumns = ceil(nAnimals/nRows);
trialSets = {'rew', 'puff', 'shock'};
for trialSet = trialSets
    trialSet = trialSet{1};
    saveName = [trialSet '_avgs_Taus'];
    ensureFigure(saveName, 1);
    for counter = 1:nAnimals
        subplot(nRows, nColumns, counter); hold on;
        plot(us_pooled.(trialSet).avgDataX{counter, 1},  us_pooled.(trialSet).avgData{counter, 1}, '--b');
        plot(us_pooled.(trialSet).avgDataX{counter, 2},  us_pooled.(trialSet).avgData{counter, 2}, '--g');
        try
            plot(us_pooled.(trialSet).tauDataX{counter, 1},  us_pooled.(trialSet).tauData{counter, 1}, '-b');
            plot(us_pooled.(trialSet).tauDataX{counter, 2},  us_pooled.(trialSet).tauData{counter, 2}, '-g');    
            t = textBox({['\color{blue}' sprintf('Ch1=%.3g',us_pooled.(trialSet).tauOff(counter,1))]; ['\color{green}' sprintf('Ch2=%.3g',us_pooled.(trialSet).tauOff(counter,2))]});             
            set(t, 'Interpreter', 'tex');
        end
        plot(us_pooled.(trialSet).latency(counter, 1), us_pooled.(trialSet).latencyY(counter, 1), '*b');
        plot(us_pooled.(trialSet).latency(counter, 2), us_pooled.(trialSet).latencyY(counter, 2), '*g');
    end
end

% saveName = 'rew_avgs_Taus';
% ensureFigure(saveName, 1);
% for counter = 1:nAnimals
%     subplot(nRows, nColumns, counter);
%     plot(xData, [gAvg.phDelay.uncuedPuff.data(counter * 2 - 1:counter*2,:) gAvg.phUs.uncuedPuff.data(counter * 2 - 1:counter*2,:)]');
% end
% 
% saveName = 'uncuedUs_avgs_shock';
% ensureFigure(saveName, 1);
% for counter = 1:nAnimals
%     subplot(nRows, nColumns, counter);
%     plot(xData, [gAvg.phDelay.uncuedShock.data(counter * 2 - 1:counter*2,:) gAvg.phUs.uncuedShock.data(counter * 2 - 1:counter*2,:)]');
% end



%% make cumulative histograms for latency and jitter for reward and
% punishment

lc = mycolors('chat');

jitter.reward = cum(us_pooled.rew.jitter_std(us_pooled.rew.avg_delta > minDelta));
jitter.puff = cum(us_pooled.puff.jitter_std(us_pooled.puff.avg_delta> minDelta));
latency.reward = cum(us_pooled.rew.latency(us_pooled.rew.avg_delta > minDelta));
latency.puff = cum(us_pooled.puff.latency(us_pooled.puff.avg_delta > minDelta));
avg.reward = cum(us_pooled.rew.avg_delta(us_pooled.rew.avg_delta > minDelta));
avg.puff= cum(us_pooled.puff.avg_delta(us_pooled.puff.avg_delta > minDelta));
tauOff.reward = cum(us_pooled.rew.tauOff(us_pooled.rew.avg_delta > minDelta));
tauOff.puff = cum(us_pooled.puff.tauOff(us_pooled.rew.avg_delta > minDelta));
figSize = [4 1];
saveName = 'cumHist_reward';
ensureFigure(saveName, 1);
subplot(1,4,1); 
plot(latency.reward.sorted * 1000, latency.reward.index, 'Color', lc); xlabel('latency (ms)'); ylabel('fraction'); set(gca, 'XLim', [0 max(get(gca, 'XLim'))]);
title('reward');
subplot(1,4,2); 
plot(jitter.reward.sorted * 1000, jitter.reward.index, 'Color', lc); xlabel('jitter (ms)'); set(gca, 'XLim', [0 max(get(gca, 'XLim'))], 'YTickLabel', {});
subplot(1,4,3); 
plot(avg.reward.sorted, avg.reward.index, 'Color', lc); xlabel('amplitude (ZS)'); set(gca, 'XLim', [0 max(get(gca, 'XLim'))], 'YTickLabel', {});
subplot(1,4,4); 
plot(tauOff.reward.sorted * 1000, tauOff.reward.index, 'Color', lc); xlabel('Tau (ms)'); set(gca, 'XLim', [0 max(get(gca, 'XLim'))], 'YTickLabel', {});
formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(figsavepath, saveName), '-eps');
end

saveName = 'cumHist_puff';
ensureFigure(saveName, 1);
subplot(1,4,1); 
plot(latency.puff.sorted * 1000, latency.puff.index, 'Color', lc); xlabel('latency (ms)'); ylabel('fraction'); set(gca, 'XLim', [0 max(get(gca, 'XLim'))]);
title('air puff');
subplot(1,4,2); 
plot(jitter.puff.sorted * 1000, jitter.puff.index, 'Color', lc); xlabel('jitter (ms)'); set(gca, 'XLim', [0 max(get(gca, 'XLim'))], 'YTickLabel', {});
subplot(1,4,3); 
plot(avg.puff.sorted, avg.puff.index, 'Color', lc); xlabel('amplitude (ZS)'); set(gca, 'XLim', [0 max(get(gca, 'XLim'))], 'YTickLabel', {});
subplot(1,4,4); 
plot(tauOff.puff.sorted * 1000, tauOff.puff.index, 'Color', lc); xlabel('Tau (ms)'); set(gca, 'XLim', [0 max(get(gca, 'XLim'))], 'YTickLabel', {}, 'YTickLabel', {});

formatFigurePublish('size', figSize);
if saveOn 
    export_fig(fullfile(figsavepath, saveName), '-eps');
end
%%
% title('appetitive');
% subplot(2,2,2); plot(nanmean([gAvgNorm.phCue.cuedReward.data gAvgNorm.phUs.cuedReward.data])); hold on; plot(nanmean([gAvgNorm.phCue.uncuedReward.data gAvgNorm.phUs.uncuedReward.data]))
% title('appetitive norm.');
% subplot(2,2,3); plot(nanmean([gAvgNorm.phCue.cuedPuff.data gAvgNorm.phUs.cuedPuff.data])); hold on; plot(nanmean([gAvgNorm.phCue.uncuedPuff.data gAvgNorm.phUs.uncuedPuff.data]))
% title('normPunish');
% subplot(2,2,4); plot(nanmean([gAvgNorm.phCue.cuedShock.data gAvgNorm.phUs.cuedShock.data])); hold on; plot(nanmean([gAvgNorm.phCue.uncuedShock.data gAvgNorm.phUs.uncuedShock.data]))
% title('normShock');


%% 
% boundedline(xData(:), [nanmean(pupData(trialsByType{3}, :)); nanmean(pupData(trialsByType{4} & TE.BlockNumber == 4, :)); nanmean(pupData(trialsByType{6}, :))]',...
%     permute([nanstd(pupData(trialsByType{3}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{3},:)), 1));...
%     nanstd(pupData(trialsByType{4} & TE.BlockNumber == 4, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{4} & TE.BlockNumber == 4,:)), 1));...
%     nanstd(pupData(trialsByType{6}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{6},:)), 1))], [2 3 1]),...
%     'cmap', [1 0 0; 0 0 0; 1 0 1]);    