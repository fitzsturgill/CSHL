savePath = 'Z:\SummaryAnalyses\ChAT_PE_Manuscript\cuedOutcome\figure\';

%%

% compile anticipatory licks and cue photometry responses for first n (50)
% trials compared to steady state

initial = {...
    'Z:\FitzRig1\Data\ChAT_26\SO_RewardPunish_odor\Session Data\', 'ChAT_26_SO_RewardPunish_odor_Jul01_2016_Session1.mat';...
    'Z:\FitzRig2\Data\ChAT_32\CuedOutcome_odor_complete\Session Data\', 'ChAT_32_CuedOutcome_odor_complete_Aug14_2016_Session1.mat';...
    'Z:\FitzRig1\Data\ChAT_34\CuedOutcome_odor_complete\Session Data\', 'ChAT_34_CuedOutcome_odor_complete_Sep23_2016_Session1.mat';...
    'Z:\FitzRig1\Data\ChAT_35\CuedOutcome_odor_complete\Session Data\', 'ChAT_35_CuedOutcome_odor_complete_Sep22_2016_Session1.mat';...
    'Z:\FitzRig2\Data\ChAT_37\CuedOutcome_odor_complete\Session Data\', 'ChAT_37_CuedOutcome_odor_complete_Sep22_2016_Session1.mat';...
    'Z:\FitzRig2\Data\ChAT_39\CuedOutcome_odor_complete\Session Data\', 'ChAT_39_CuedOutcome_odor_complete_Sep22_2016_Session1.mat';...
    'Z:\FitzRig1\Data\ChAT_42\CuedOutcome_odor_complete\Session Data\', 'ChAT_42_CuedOutcome_odor_complete_Sep22_2016_Session1.mat'...
    };

DB = dbLoadExperiment('cuedOutcome');
DBSO = dbLoadExperiment('SO_RewardPunish_odor');

firstTrials = 25; % because you are just using the first X trials, don't worry about truncating this first session for now, you are explicitly taking the beginning of the session

data = struct(...
    'initial', [],...
    'others', []...
    );
pooled = struct(...
    'csLicks', data,...
    'csFluor', data...
    );
    % Note! for SO_RewardPunish_odor, trial types 1-2 comprise the
    % appetitive odor, for cuedOutcome_Odor_Complete, trial types 1-3 do
    % the same
    window = [-3 6];
    [earlyPh, earlyLicks] = deal(zeros(40, diff(window) * 20, 7));
    
for counter = 1:length(DB.animals)
    sessions = bpLoadSessions([], initial{counter,2}, initial{counter,1});
    if counter == 1 % this is a SO_rewardpunish_odor mouse
        dbLoadAnimal(DBSO, DB.animals{counter});
        TE.phPeak_cs = bpCalcPeak_dFF(TE.Photometry, 1, [0 2], TE.Cue, 'method', 'mean', 'phField', 'ZS');        
        
        TE_temp = makeTE_SO_RewardPunish_Odor_v2(sessions);
        TE_temp.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', 1, 'baseline', [0 4]);        
        TE_temp.csLicks = countEventFromTE(TE_temp, 'Port1In', [0 2], TE_temp.Cue);
        TE_temp.usLicks = countEventFromTE(TE_temp, 'Port1In', [0 1], max([cellfun(@(x) x(1), TE_temp.Reward) cellfun(@(x) x(1), TE_temp.Punish) cellfun(@(x) x(1), TE_temp.Omit)], 2));
        TE_temp.phPeak_cs = bpCalcPeak_dFF(TE_temp.Photometry, 1, [0 2], TE_temp.Cue, 'method', 'mean', 'phField', 'ZS');
        keepers = find(ismember(TE_temp.trialType, [1 2]));
%       phAverageFromTE(TE_temp, keepers(1:min(firstTrials, length(keepers))), 1, 'zeroTimes', TE_temp.Cue, 'window', window, 'FluorDataField', 'ZS');
        [data, xData] = phAlignedWindow(TE_temp, keepers(1:40), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE_temp.Cue, 'window', window);
        earlyPh(:,:,counter) = data;
        eventTS = bpEventToTimeSeries(TE_temp, 'Port1In', 'zeroField', 'Cue', 'duration', 10);
        data = alignedDataWindow(eventTS.data, keepers(1:40), 'zeroTimes', TE_temp.Cue, 'window', window, 'Fs', 20, 'startTimes', TE_temp.Photometry.startTime);
        earlyLicks(:,:,counter) = data;
        
        pooled.csLicks.initial{counter} = TE_temp.csLicks.rate(keepers(1:min(firstTrials, length(keepers))));
        pooled.csFluor.initial{counter} = TE_temp.phPeak_cs.data(keepers(1:min(firstTrials, length(keepers))));
        pooled.csFluor.others{counter} = TE.phPeak_cs.data(rewardOdorTrials);
        pooled.csLicks.others{counter} = TE.csLicks.rate(rewardOdorTrials);
    else % these are cuedOutcome_odor_complete mice
        dbLoadAnimal(DB, DB.animals{counter});
        TE_temp = makeTE_CuedOutcome_Odor_Complete(sessions);
        TE_temp.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', 1, 'baseline', [0 4]);
        TE_temp.csLicks = countEventFromTE(TE_temp, 'Port1In', [0 2], TE_temp.Cue);
        TE_temp.usLicks = countEventFromTE(TE_temp, 'Port1In', [0 1], TE_temp.Us);
        TE_temp.phPeak_cs = bpCalcPeak_dFF(TE_temp.Photometry, 1, [0 2], TE_temp.Cue, 'method', 'mean', 'phField', 'ZS');        
        keepers = find(ismember(TE_temp.trialType, [1 2 3]));
        pooled.csLicks.initial{counter} = TE_temp.csLicks.rate(keepers(1:min(firstTrials, length(keepers))));
        pooled.csFluor.initial{counter} = TE_temp.phPeak_cs.data(keepers(1:min(firstTrials, length(keepers))));
        pooled.csFluor.others{counter} = TE.phPeak_cs.data(highValueTrials);
        pooled.csLicks.others{counter} = TE.csLicks.rate(highValueTrials);
        [data, ~] = phAlignedWindow(TE_temp, keepers(1:40), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE_temp.Cue, 'window', window);
        earlyPh(:,:,counter) = data;
        
        eventTS = bpEventToTimeSeries(TE_temp, 'Port1In', 'zeroField', 'Cue', 'duration', 10);
        data = alignedDataWindow(eventTS.data, keepers(1:40), 'zeroTimes', TE_temp.Cue, 'window', window, 'Fs', 20, 'startTimes', TE_temp.Photometry.startTime);
        earlyLicks(:,:,counter) = data;
        
%         avgData = phAverageFromTE(TE, highValueTrials, 1, 'zeroTimes', TE.Cue, 'window', window, 'FluorDataField', 'ZS');
%         lateCue(counter,:) = avgData.Avg;
    end
end

%%
ensureFigure('test', 1);
meanData = [cellfun(@nanmean, pooled.csFluor.initial); cellfun(@nanmean, pooled.csFluor.others)];
plot(meanData)

%%
ensureFigure('test', 1);
for counter = 1:7
    subplot(3,3,counter);
    plot(mean(earlyLicks(:,:,counter)));
end
%% averages for trials 1-10,11-20,21-30,31-40 from first sessions
earlyAvgs = permute(earlyPh,[2 3 1]);
earlyAvgs = reshape(earlyAvgs, 180, 10*7, 4);
earlyAvgs = squeeze(median(earlyAvgs, 2));

% get rid of outliers
outliers = earlyLicks > 5;
earlyLicks(outliers) = 5;
earlyLickAvgs = permute(earlyLicks, [2 3 1]);
earlyLickAvgs = reshape(earlyLickAvgs, 180, 10*7, 4);
earlyLickAvgs = squeeze(mean(earlyLickAvgs, 2));

temp = cool;
cmap = temp(round(linspace(10,64, 4)), :);


saveName = 'grandAvgs_initial_novelty';
ensureFigure(saveName, 1);
axes; hold on;
for counter = 1:4
    plot(xData, smooth(earlyAvgs(:,counter)), 'Color', cmap(counter, :));
end

addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);

xlabel('time from odor (s)'); ylabel('Fluor. (\sigma-bl.)');
set(gca, 'XLim', [-2 6], 'YTick', [0 1 2]);
formatFigurePublish('size', [2 1]);    

if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
end    

saveName = 'grandAvgs_initial_novelty_licks';
ensureFigure(saveName, 1);
axes; hold on;
for counter = 1:4
    plot(xData, smooth(earlyLickAvgs(:,counter)), 'Color', cmap(counter, :));
end

addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);

xlabel('time from odor (s)'); ylabel('licks/s');
set(gca, 'XLim', [-2 6], 'YTick', [0 1]);
formatFigurePublish('size', [2 1]);    

if saveOn
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
end    


%% collect and average data, spaced by every n trials to show lack of novelty coding
nTrials = 10;

% figure out max numbers of trials, round up to be divisible by nTrials
maxPoints = max(cellfun(@length, pooled.csFluor.initial) + cellfun(@length, pooled.csFluor.others));
maxPoints = ceil(maxPoints / nTrials) * nTrials;
[fluorData, lickData] = deal(NaN(maxPoints/nTrials, 7));

for counter = 1:7
    data = [pooled.csFluor.initial{counter}; pooled.csFluor.others{counter}];
    nPoints = length(data);
    data = [data; NaN(maxPoints - nPoints,1)];
    data = reshape(data, nTrials, maxPoints/nTrials);
    fluorData(:,counter) = nanmean(data)';
    
    data = [pooled.csLicks.initial{counter}; pooled.csLicks.others{counter}];
    data = [data; NaN(maxPoints - nPoints,1)];
    data = reshape(data, nTrials, maxPoints/nTrials);
    lickData(:,counter) = nanmean(data)';
end

nShow = 10;
figSize = [2 1];

saveName = 'novelty_bar_fluor';
ensureFigure(saveName, 1);
axes;

xData = 1:nTrials:(nShow*nTrials);
yData = nanmean(fluorData(1:nShow,:),2);
errs = nanSEM(fluorData(1:nShow,:)');

errorbar(xData(:), yData(:), errs(:), 'k');
formatFigurePublish('size', figSize);
set(gca, 'YTick', [0 1]);
xlabel('odor trials');
ylabel('Fluor');

if saveOn 
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));    
end

saveName = 'novelty_bar_lick';
ensureFigure(saveName, 1);
axes;

xData = 1:nTrials:(nShow*nTrials);
yData = nanmean(lickData(1:nShow,:),2);
errs = nanSEM(lickData(1:nShow,:)');

errorbar(xData(:), yData(:), errs(:), 'k');
formatFigurePublish('size', figSize);
xlabel('odor trials');
ylabel('licks/s');

if saveOn 
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));    
end

%% example lick and photometry raster for lack of novelty coding
% #6 corresponds to ChAT_39
    
window = [-4 6];
sessions = bpLoadSessions([], initial{6,2}, initial{6,1});
    

TE_temp = makeTE_CuedOutcome_Odor_Complete(sessions);
TE_temp.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', 1, 'baseline', [0 4]);
keepers = find(ismember(TE_temp.trialType, [1 2 3]));
%%
figSize = [1.2 1.6];
climfactor = 2.5;

saveName = 'novelty_raster_lick';
ensureFigure(saveName, 1);
[~, lh] = eventRasterFromTE(TE_temp, keepers, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
set(gca, 'XLim', window);
set(lh, 'LineWidth', 0.3, 'Color', [0 0 0]);
formatFigurePublish('size', figSize);
xlabel('time from odor (s)');
ylabel('trial #');

if saveOn 
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));    
end

saveName = 'novelty_raster_fluor';
ensureFigure(saveName, 1);
phRasterFromTE(TE_temp, keepers, 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', 'ZS', 'PhotometryField', 'Photometry', 'zeroTimes', TE_temp.Cue, 'window', window); % 'CLimFactor', CLimFactor,
set(gca, 'XLim', window);

formatFigurePublish('size', figSize);
xlabel('time from odor (s)');
ylabel('trial #');

if saveOn 
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));    
end

    
    
    
        

































%%
% ensureFigure('test2', 1);
% for counter = 1:7
%     subplot(4,2,counter); hold on;
%     nums = [length(pooled.csFluor.initial{counter}) length(pooled.csFluor.others{counter})];
%     plot(1:nums(1), pooled.csFluor.initial{counter}, 'g');
%     plot(nums(1)+1:nums(1) + nums(2), pooled.csFluor.others{counter}, 'k');
% end
% 
% ensureFigure('csLicks', 1);
% for counter = 1:7
%     subplot(4,2,counter); hold on;
%     nums = [length(pooled.csLicks.initial{counter}) length(pooled.csLicks.others{counter})];
%     plot(1:nums(1), pooled.csLicks.initial{counter}, 'g');
%     plot(nums(1)+1:nums(1) + nums(2), pooled.csLicks.others{counter}, 'k');
% end
