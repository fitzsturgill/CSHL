%



DB = dbLoadExperiment('FrankenLNL_RewardPunish');
figPath = fullfile(DB.path, 'figure');
ensureDirectory(figPath);
smoothWindow = 1;
saveOn = 1;
figSize = [2 1];

%% examples for ACh_7 and ACh_15, showing similar signals, but different degrees of noise correlations

photometryField = 'Photometry';
fdField = 'ZS';
boundField = 'SEM';
window = [-6.5 4];
linecolors = [mycolors('reward'); mycolors('reward_cued')];
lineWidth = 0.25;

animal = 'ACh_7';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups
saveName = sprintf('PE_BLA_%s_avgs', animal);  
h=ensureFigure(saveName, 1); 
sessionIndexList = 5; % just the 1 session because early sessions surprise modulation hasn't developed whereas later sessions I shorten the delay (surprise modulation is consistent) but it messes up the graph having 2 different delays


avgData1 = phAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 1,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window);
avgData2 = phAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 2,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window);
rch = [max(range(avgData1.Avg')) max(range(avgData2.Avg'))];
offset = rch(1) + rch(2) * 0.25;
xData = avgData1.xData(1,:);

subplot(1,2,2); hold on;
[thisHl, thisHp] = boundedline(xData, avgData1.Avg, permute(avgData1.(boundField), [2 3 1]), gca, 'cmap', linecolors, 'alpha', 'nan', 'gap');
set(thisHl, 'LineWidth', lineWidth);
[thisHl, thisHp] = boundedline(xData, avgData2.Avg + offset, permute(avgData2.(boundField), [2 3 1]), gca, 'cmap', linecolors, 'alpha', 'nan', 'gap');
set(thisHl, 'LineWidth', lineWidth);
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
set(gca, 'XLim', window);
xlabel('Time (s)');

animal = 'ACh_15';
sessionIndexList = [6 7];
success = dbLoadAnimal(DB, animal); % load TE and trial lookups
avgData1 = phAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 1,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window);
avgData2 = phAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 2,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window);
avgData1.Avg = avgData1.Avg * 1;
avgData1.(boundField) = avgData1.(boundField) * 1;
rch = [max(range(avgData1.Avg')) max(range(avgData2.Avg'))];
offset = rch(1) + rch(2) * 0.25;
xData = avgData1.xData(1,:);

subplot(1,2,1); hold on;
[thisHl, thisHp] = boundedline(xData, avgData1.Avg, permute(avgData1.(boundField), [2 3 1]), gca, 'cmap', linecolors, 'alpha', 'nan', 'gap');
set(thisHl, 'LineWidth', lineWidth);
[thisHl, thisHp] = boundedline(xData, avgData2.Avg + offset, permute(avgData2.(boundField), [2 3 1]), gca, 'cmap', linecolors, 'alpha', 'nan', 'gap');
set(thisHl, 'LineWidth', lineWidth);
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
set(gca, 'XLim', window);
xlabel('Time (s)');
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');

formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(figPath, [saveName '.pdf']));
    saveas(gcf, fullfile(figPath, [saveName '.fig']));
    saveas(gcf, fullfile(figPath, [saveName '.jpg'])); 
end

%% plot grand averages for cued, uncued and omitted reward, use only block 2, introductory block, trace = 1s
% unlike RankenLNL_RewardPunish_poolAnimals, focus on more restricted
% training phase to ensure a constant trace delay

DB = dbLoadExperiment('FrankenLNL_RewardPunish');
saveOn = 1;
minRewardLickRate = 2; % at least n Hz licking (during us window for reward trials)
minCueLickRate = 1; % at least n Hz licking (during cs window for reward cued reward trials)
savepath = fullfile(DB.path, ['pooled' filesep]);
ensureDirectory(savepath);
figsavepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(figsavepath);

% Goal 1: Grand Averages
% data and associated descriptors for each data type, trial set combination
s3 = struct(...
    'data', [],...
    'animal', [],...
    'ch', []...
    );

% initialize different trial setsa
s2 = struct(...
    'cuedReward', s3,...
    'uncuedReward', s3,...
    'omitReward', s3...
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
    animal = DB.animals{counter};
%     if strcmp(animal, 'ACh_3')
%         continue;
%     end
    dbLoadAnimal(DB, animal);

    % for each animal, assemble trial sets, do this manually to allow tweaking of trial makeup (exclude conditions such as late licking, no licking, select certain block numbers, etc.)
    %first column matches fields in gAvg, second column contains trial vector
    trialSets = {};
    
    % cuedReward
    trialSets{end+1, 1} = 'cuedReward';
    trialSets{end,2} = rewardTrials & Odor2Valve1Trials & (TE.licks_us.rate > minRewardLickRate) & (TE.licks_cs.rate > minCueLickRate) & ismember(TE.BlockNumber, [2])...
        & (round(cellfun(@(x) diff(x), TE.Trace2)) == 1); % exclude shock or extinction days
    trialSets{end,3} = 2;
    % uncuedReward
    trialSets{end+1, 1} = 'uncuedReward';
    trialSets{end,2} = uncuedReward & (TE.licks_us.rate > minRewardLickRate) & ismember(TE.BlockNumber, [2])...
        & (round(cellfun(@(x) diff(x), TE.Trace2)) == 1);    
    trialSets{end,3} = 2;    
    % omitReward
    trialSets{end+1, 1} = 'omitReward';
    trialSets{end,2} = neutralTrials & Odor2Valve1Trials & ismember(TE.BlockNumber, [2])...
        & (round(cellfun(@(x) diff(x), TE.Trace2)) == 1);        
    trialSets{end,3} = 2;    
   
    
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

% compute averages, SEM, also normalize grand Average data to uncued reward us amplitude
gAvgNorm = gAvg;
% photometry reward Us size per fiber
normVectorPh = max(gAvg.phUs.uncuedReward.data(:,1:bpX2pnt(1, 20)), [], 2);
% lick reward Us size per fiber
% normVectorLick = nanmean(gAvg.lickUs.uncuedReward.data(:,1:bpX2pnt(1, 20)), 2);
for tcounter = 1:size(trialSets, 1)
    label = trialSets{tcounter, 1};
    gAvg.phCue.(label).Avg = nanmean(gAvg.phCue.(label).data);
    gAvg.phCue.(label).SEM = nanstd(gAvg.phCue.(label).data, 0, 1) ./ sqrt(sum(isfinite(gAvg.phCue.(label).data), 1));
    gAvg.phCue.(label).xData = (0:(size(gAvg.phCue.(label).data, 2) - 1)) * 1/20 - 7;
    gAvg.phDelay.(label).Avg = nanmean(gAvg.phDelay.(label).data);
    gAvg.phDelay.(label).SEM = nanstd(gAvg.phDelay.(label).data, 0, 1) ./ sqrt(sum(isfinite(gAvg.phDelay.(label).data), 1));
    gAvg.phDelay.(label).xData = (0:(size(gAvg.phDelay.(label).data, 2) - 1)) * 1/20 - 2;    
    gAvg.ph.(label).Avg = nanmean(gAvg.ph.(label).data);
    gAvg.ph.(label).SEM = nanstd(gAvg.ph.(label).data, 0, 1) ./ sqrt(sum(isfinite(gAvg.ph.(label).data), 1));
    gAvg.ph.(label).xData = (0:(size(gAvg.ph.(label).data, 2) - 1)) * 1/20 - (size(gAvg.ph.(label).data, 2)/20 - 4);
    gAvg.phUs.(label).Avg = nanmean(gAvg.phUs.(label).data);
    gAvg.phUs.(label).SEM = nanstd(gAvg.phUs.(label).data, 0, 1) ./ sqrt(sum(isfinite(gAvg.phUs.(label).data), 1));    
    gAvg.phUs.(label).xData = (0:(size(gAvg.phUs.(label).data, 2) - 1)) * 1/20;
    
    gAvgNorm.phCue.(label).data = gAvg.phCue.(label).data ./ normVectorPh;
    gAvgNorm.phDelay.(label).data = gAvg.phDelay.(label).data ./ normVectorPh;
    gAvgNorm.ph.(label).data = gAvg.ph.(label).data ./ normVectorPh;
    gAvgNorm.phUs.(label).data = gAvg.phUs.(label).data ./ normVectorPh;
    
    gAvgNorm.phCue.(label).Avg = nanmean(gAvgNorm.phCue.(label).data);
    gAvgNorm.phCue.(label).SEM = nanstd(gAvgNorm.phCue.(label).data, 0, 1) ./ sqrt(sum(isfinite(gAvgNorm.phCue.(label).data), 1));
    gAvgNorm.phCue.(label).xData = (0:(size(gAvgNorm.phCue.(label).data, 2) - 1)) * 1/20 - 7;
    gAvgNorm.phDelay.(label).Avg = nanmean(gAvg.phDelay.(label).data);
    gAvgNorm.phDelay.(label).SEM = nanstd(gAvg.phDelay.(label).data, 0, 1) ./ sqrt(sum(isfinite(gAvg.phDelay.(label).data), 1));
    gAvgNorm.phDelay.(label).xData = (0:(size(gAvg.phDelay.(label).data, 2) - 1)) * 1/20 - 2;        
    gAvgNorm.ph.(label).Avg = nanmean(gAvgNorm.ph.(label).data);
    gAvgNorm.ph.(label).SEM = nanstd(gAvgNorm.ph.(label).data, 0, 1) ./ sqrt(sum(isfinite(gAvgNorm.ph.(label).data), 1));
    gAvgNorm.ph.(label).xData = (0:(size(gAvgNorm.ph.(label).data, 2) - 1)) * 1/20 - (size(gAvgNorm.ph.(label).data, 2)/20 - 4);
    gAvgNorm.phUs.(label).Avg = nanmean(gAvgNorm.phUs.(label).data);
    gAvgNorm.phUs.(label).SEM = nanstd(gAvgNorm.phUs.(label).data, 0, 1) ./ sqrt(sum(isfinite(gAvgNorm.phUs.(label).data), 1));    
    gAvgNorm.phUs.(label).xData = (0:(size(gAvgNorm.phUs.(label).data, 2) - 1)) * 1/20;
    
%     gAvgNorm.lickCue(label).data = gAvgNorm.lickCue(label).data ./ normVectorLick;
%     gAvgNorm.lick(label).data = gAvgNorm.lick(label).data ./ normVectorLick;
%     gAvgNorm.lickUs(label).data = gAvgNorm.lickUs(label).data ./ normVectorLick;
end




%% plot grand averages

figSize = [1.7 0.9];
saveName = 'grandAverage_appetitive_simple';
ensureFigure(saveName, 1); 
linecolors = [mycolors('reward'); 0 0 0; mycolors('reward_cued')];     
window = [-4 4];
axes; hold on;
% xData = [gAvg.phCue.cuedReward.xData gAvg.phUs.cuedReward.xData;...
%     gAvg.phCue.omitReward.xData gAvg.phUs.omitReward.xData;...
%     gAvg.phCue.uncuedReward.xData gAvg.phUs.uncuedReward.xData]';
% yData = [gAvg.phCue.cuedReward.Avg gAvg.phUs.cuedReward.Avg;...
%     gAvg.phCue.omitReward.Avg gAvg.phUs.omitReward.Avg;
%     gAvg.phCue.uncuedReward.Avg gAvg.phUs.uncuedReward.Avg]';
% bData = permute([gAvg.phCue.cuedReward.SEM gAvg.phUs.cuedReward.SEM;...
%     gAvg.phCue.omitReward.SEM gAvg.phUs.omitReward.SEM;...
%     gAvg.phCue.uncuedReward.SEM gAvg.phUs.uncuedReward.SEM], [2 3 1]);

xData = [gAvg.ph.cuedReward.xData;...
    gAvg.ph.omitReward.xData;...
    gAvg.ph.uncuedReward.xData]';
yData = [gAvg.ph.cuedReward.Avg;...
    gAvg.ph.omitReward.Avg;
    gAvg.ph.uncuedReward.Avg]';
bData = permute([gAvg.ph.cuedReward.SEM;...
    gAvg.ph.omitReward.SEM;...
    gAvg.ph.uncuedReward.SEM], [2 3 1]);

[hl, hp] = boundedline(xData, yData, bData, 'cmap', linecolors);
% set(gca, 'XLim', window);
set(gca, 'YLim', [-0.8 2.5]);
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% legend(hvl, {'cued', 'omit', 'uncued'}, 'Location', 'best'); legend('boxoff');
ylabel('F(\fontsize{10}\sigma\fontsize{7}-baseline)');  set(gca, 'XLim', window);
xlabel('Time from reinforcement (s)');

formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(figsavepath, [saveName '.pdf']));
    export_fig(fullfile(figsavepath, saveName), '-eps');
end






