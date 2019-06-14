%{
developed as subset to reversals_noPunish_poolReversals, focuses on trying
to leverage statistical power of paired recordings of ACh. and Dop. neurons

%}
DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, 'pooled', filesep);
smoothWindow = 5; % smooth data to make tt50 less sensitive to noise
saveOn = 1;
baselineTrials = 20;

%%
exp.value = {'DC_44'  'DC_46'  'DC_47' 'DC_53'  'DC_54'  'DC_56'}; % exclude DC_51
exp.valence = {'DC_17'  'DC_20'  'DC_35'  'DC_36'  'DC_37'  'DC_40'};
exp.all = [exp.value exp.valence];
expType = 'all';


%%
compile_reversal_data;

%% HACK ALERT- TO INCLUDE PUNISH REVERSALS TEMPORARILY (or not so temporarily)
% HACK ALERT- TO INCLUDE PUNISH REVERSALS TEMPORARILY
% HACK ALERT- TO INCLUDE PUNISH REVERSALS TEMPORARILY
% HACK ALERT- TO INCLUDE PUNISH REVERSALS TEMPORARILY
% HACK ALERT- TO INCLUDE PUNISH REVERSALS TEMPORARILY

rocThresh = 0.5;
trialsToCriterion = NaN(nReversals, 1);
% looping is just easier
for counter = 1:nReversals    
    thisRev = AR.csPlus.csLicksROC.after(counter, :);
    thisRev = thisRev > rocThresh;
    nt = find(thisRev, 1);
    if ~isempty(nt) && any(isfinite(AR.csPlus.csLicksROC.after(counter, :)))
        trialsToCriterion(counter) = nt;
    else
        trialsToCriterion(counter) = Inf;  % HACK
    end
end

%% filter reversals according to quality
goodReversals = ...
    ~isnan(trialsToCriterion) &...
    auROC.phPeakMean_cs_ch1.acq > 0 &...
    auROC.phPeakMean_cs_ch2.acq > 0;

% goodReversals = ...
%     auROC.phPeakMean_cs_ch1.acq > 0 &...
%     auROC.phPeakMean_cs_ch2.acq > 0 &...
%     auROC.licks_cs.after > 0;...

sortVariable = trialsToCriterion;
% sortVariable = cp.csPlus.phPeakMean_cs_ch1.index;
sortVariable(~goodReversals) = NaN;

[sorted, sortOrder] = sort(sortVariable);
% sortOrder = sortOrder(~isnan(sorted));


%% take advantage of paired recordings, subtract dopamine cue response from Ach. cue response, etc.
newCsPlus.phPeakMean_cs_AchMinusDop = newCsPlus.phPeakMean_cs_ch1 - newCsPlus.phPeakMean_cs_ch2;
newCsMinus.phPeakMean_cs_AchMinusDop = newCsMinus.phPeakMean_cs_ch1 - newCsMinus.phPeakMean_cs_ch2;
alwaysCsPlus.phPeakMean_cs_AchMinusDop = alwaysCsPlus.phPeakMean_cs_ch1 - alwaysCsPlus.phPeakMean_cs_ch2;
odor3.phPeakMean_cs_AchMinusDop = odor3.phPeakMean_cs_ch1 - odor3.phPeakMean_cs_ch2;

%% find trials to 50%
fraction = 0.5;
compFields = {'csPlus', 'csMinus';...  % first row is after reversal
              'csMinus', 'csPlus'...   % second row is before reversal
              };
fitFields = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'};
outputFields = {'bl', 'max', 'tt50'};
tt50 = struct();

nReversals = sum(goodReversals);

for counter = 1:length(compFields)
    for counter2 = 1:length(fitFields)
        for counter3 = 1:length(outputFields)            
            tt50.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = NaN(nReversals, 1);
        end
    end
end
tt50.fraction = fraction;
tt50.baselineTrials = baselineTrials;


for compCounter = 1:size(compFields, 2)    
    for fieldCounter = 1:length(fitFields)
        bl = nanmean(AR.(compFields{2, compCounter}).(fitFields{fieldCounter}).before(goodReversals, end - baselineTrials + 1:end), 2);
        tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).bl = bl;
        tt50Data = AR.(compFields{1, compCounter}).(fitFields{fieldCounter}).after(goodReversals, :);
        for counter = 1:size(tt50Data, 1)        
            thisData = tt50Data(counter, :);
            if ~any(isfinite(thisData))
                continue
            end
            switch compFields{1, compCounter}
                case 'csPlus'
                    top = percentile(thisData, 0.9);
                    thresh = (bl(counter) + (top - bl(counter)) * fraction);
                    latency = find(thisData > thresh, 1);
                    if isempty(latency)
                        continue
                    end
                    tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).tt50(counter) = latency;
                    tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).thresh(counter) = thresh;
                case 'csMinus'
                    bottom = percentile(thisData, 0.1);
                    thresh = (bl(counter) - (bl(counter) - bottom) * fraction);
                    latency = find(thisData < thresh, 1);
                    if isempty(latency)
                        continue
                    end
                    tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).tt50(counter) = latency;                    
                    tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).thresh(counter) = thresh;
            end            
        end
    end
end


% repeat, but permute data
% randomly permute across reversals (not time) to preserve average but
% scramble potential tt50s
maxTrials = size(AR.csPlus.filename.after, 2);

% generate subscript indices for permutation
row_ix = repmat(1:maxTrials, sum(goodReversals), 1);
col_ix = zeros(sum(goodReversals), maxTrials);
for counter = 1:maxTrials
    col_ix(:,counter) = randperm(sum(goodReversals))';
end

good_licks = AR.csPlus.licks_cs.after(goodReversals, :);
good_subtract = AR.csPlus.phPeakMean_cs_ch1.after(goodReversals, :) - AR.csPlus.phPeakMean_cs_ch2.after(goodReversals, :);
good_ch1 = AR.csPlus.phPeakMean_cs_ch1.after(goodReversals, :);
good_ch2 = AR.csPlus.phPeakMean_cs_ch2.after(goodReversals, :);

lin_ix = sub2ind([sum(goodReversals) maxTrials], col_ix, row_ix);
perm_subtract = smoothdata(good_subtract(lin_ix), 2, 'movmean', smoothWindow, 'omitnan');
perm_licks = smoothdata(good_licks(lin_ix), 2, 'movmean', smoothWindow, 'omitnan');
perm_ch1 = smoothdata(good_ch1(lin_ix), 2, 'movmean', smoothWindow, 'omitnan');
perm_ch2 = smoothdata(good_ch2(lin_ix), 2, 'movmean', smoothWindow, 'omitnan');

% quick plot of regular and permuted lick_cs images, and corresponding
% averages
savename = 'permutation_sanity_check';
ensureFigure(savename, 1);
clim = [0 20]; xlim = [0 40];
subplot(3,1,1);
imagesc(good_licks, clim);
set(gca, 'XLim', xlim);
title('Licks');
subplot(3,1,2);
imagesc(perm_licks, clim);
set(gca, 'XLim', xlim);
title('Licks permuted');
subplot(3,1,3); hold on; 
title('averages');
plot(nanmean(good_licks), '-k');
plot(nanmean(perm_licks), '--', 'Color', [0.8 0.8 0.8]);
legend({'intact', 'permuted'}, 'Location', 'best');
set(gca, 'XLim', xlim); xlabel('new Cs+ trials from rev.');

formatFigurePoster([3 5], '', 10);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
end    


% detect tt50 on permuted data
tt50_perm = struct();


tt50_perm.csPlus.licks_cs = NaN(sum(goodReversals), 1);
tt50_perm.fraction = fraction;
tt50_perm.baselineTrials = baselineTrials;


% just use a scalar baseline value for simplicity here
bl = AR.csPlus.licks_cs.before(goodReversals, end - baselineTrials + 1:end);
bl = bl(:);
bl = nanmean(bl);

tt50_perm.csPlus.bl = bl;
tt50Data = perm_licks;
for counter = 1:size(tt50Data, 1)        
    thisData = tt50Data(counter, :);
    if ~any(isfinite(thisData))
        continue
    end
%     switch compFields{1, compCounter}
%         case 'csPlus'
            top = percentile(thisData, 0.9);
            thresh = (bl + (top - bl) * fraction);
            latency = find(thisData > thresh, 1);
            if isempty(latency)
                continue
            end
            tt50_perm.csPlus.tt50(counter) = latency;
            tt50_perm.csPlus.thresh(counter) = thresh;
%         case 'csMinus'
%             bottom = percentile(thisData, 0.1);
%             thresh = (bl(counter) - (bl(counter) - bottom) * fraction);
%             latency = find(thisData < thresh, 1);
%             if isempty(latency)
%                 continue
%             end
%             tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).tt50(counter) = latency;                    
%             tt50.(compFields{1, compCounter}).(fitFields{fieldCounter}).thresh(counter) = thresh;
%     end            
end




%% make a bar graph or box plot of time to 50% (or other fraction)
    
savename = ['tt50_all_' expType];
ensureFigure(savename, 1); axes('FontSize', 12); hold on;
markerSize = 15;

fields = {...
    'csPlus', 'licks_cs', mycolors('licks'), 'licks';...
    'csPlus', 'phPeakMean_cs_ch1', mycolors('chat'), 'ACh.';...
    'csPlus', 'phPeakMean_cs_ch2', mycolors('dat'), 'Dop.';...
    'csMinus', 'licks_cs', mycolors('licks'), 'licks';...
    'csMinus', 'phPeakMean_cs_ch1', mycolors('chat'), 'ACh.';...
    'csMinus', 'phPeakMean_cs_ch2', mycolors('dat'), 'Dop.';...    
    };


all_tt50 = [];

for counter = 1:size(fields, 1)
    ydata = tt50.(fields{counter, 1}).(fields{counter, 2}).tt50;       
    all_tt50 = [all_tt50 ydata];
%     scatter(repmat(counter, numel(ydata), 1) + (rand(numel(ydata), 1) - 0.5)/3, ydata, markerSize, fields{counter, 3}, '.');    
%     errorbar(counter, mean(ydata), std(ydata)/sqrt(numel(ydata)), 'Color', fields{counter, 3}, 'LineWidth', 2)
end
violins = violinplot(all_tt50, fields(:, 4), 'ShowMean', true, 'ShowNotches', true');
for counter = 1:length(violins)
    violins(1).ViolinColor = fields{counter,3};
end
ylabel(sprintf('trials to %.2g%% max cue licks', fraction * 100));
% set(gca, 'XLim', [0.5 6.5], 'XTick', 1:6, 'XTickLabel', fields(:,4)', 'FontSize', 12); ylabel('trials from rev.', 'FontSize', 12);
set(gca, 'YLim', [0 20]);

formatFigurePublish('size', [3.5 2.5], 'fontSize', 12);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    export_fig(fullfile(savepath, savename), '-eps');
end






%% align reversals by tt50

trialWindow = [-20 20];
cpField = 'licks_cs';
cLimFactor = 3;
trialWindow_rev = [0 40];
markerSize = 4;
smoothWindow = 1;


% first new csPlus (acquisition)

% setup for aligned by changepoint
goodOnes = goodReversals;
zeroTrials = tt50.csPlus.(cpField).tt50;
reversalPoints = 0 - zeroTrials;
reversalPoints = reversalPoints;
[sorted, ix] = sort(reversalPoints); % THEN sort them


trials = true(sum(goodOnes), 1);
[aligned_subtract, xData] = alignedDataWindow(good_subtract, trials, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', ones(sum(goodOnes), 1));
[aligned_licks, ~] = alignedDataWindow(good_licks, trials, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', ones(sum(goodOnes), 1));
[aligned_ch1, ~] = alignedDataWindow(good_ch1, trials, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', ones(sum(goodOnes), 1));
[aligned_ch2, ~] = alignedDataWindow(good_ch2, trials, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', ones(sum(goodOnes), 1));




zeroTrials_perm = tt50_perm.csPlus.tt50;
reversalPoints_perm = 0 - zeroTrials_perm;
[sorted_perm, ix_perm] = sort(reversalPoints_perm); % THEN sort them

[aligned_subtract_perm, xData] = alignedDataWindow(perm_subtract, trials, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', ones(sum(goodOnes), 1));
[aligned_licks_perm, ~] = alignedDataWindow(perm_licks, trials, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', ones(sum(goodOnes), 1));
[aligned_ch1_perm, ~] = alignedDataWindow(perm_ch1, trials, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', ones(sum(goodOnes), 1));
[aligned_ch2_perm, ~] = alignedDataWindow(perm_ch2, trials, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', ones(sum(goodOnes), 1));


% setup for alignment by reversal
% [sorted_rev, sortOrder_rev] = sort(tt50.csPlus.licks_cs.tt50(goodOnes));
tt50_rev = zeroTrials;
tt50_rev_perm = zeroTrials_perm;

% first images, sorted by reversal point.
savename = ['tt50_aligned_newCsPlus_images' '_' expType];
ensureFigure(savename, 1);



subplot(4,4,1); hold on;
title('Reversal Aligned');
ylabel('Cue licks');
cData = good_licks;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc(cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(tt50_rev, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);

subplot(4,4,2); hold on;
title('TT50 Aligned');
cData = aligned_licks(ix, :);
imagesc('XData', xData,  'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4);
set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,3); hold on;
title('Rev. aligned, permuted');
cData = perm_licks;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc(cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(tt50_rev_perm, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);


subplot(4,4,4); hold on;
title('TT50 aligned, permuted');
cData = aligned_licks_perm(ix_perm, :);
imagesc('XData', xData,  'CData', cData, clim);
plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,5); hold on;
ylabel('ACh.', 'Color', mycolors('ChAT'));
cData = good_ch1;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc(cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(tt50_rev, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

subplot(4,4,6); hold on;
cData = aligned_ch1(ix, :);
imagesc('XData', xData,  'CData',  cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,7); hold on;
cData = perm_ch1;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc(cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(tt50_rev_perm, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);

subplot(4,4,8); hold on;
cData = aligned_ch1_perm(ix_perm, :);
imagesc('XData', xData,  'CData',  cData, clim);
plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,9); hold on;
ylabel('Dop.', 'Color', mycolors('DAT'));
cData = good_ch2;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc(cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(tt50_rev, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

subplot(4,4,10); hold on;
cData = aligned_ch2(ix, :);
imagesc('XData', xData,  'CData',  cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,11); hold on;
cData = perm_ch2;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc(cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(tt50_rev_perm, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);

subplot(4,4,12); hold on;
cData = aligned_ch2_perm(ix_perm, :);
imagesc('XData', xData,  'CData',  cData, clim);
plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,13); hold on;
ylabel('Subtract', 'Color', [0 1 0]);
cData = good_subtract;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc(cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(tt50_rev, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);
xlabel('Trials from reversal');

subplot(4,4,14); hold on;
cData = aligned_subtract(ix, :);
imagesc('XData', xData,  'CData',  cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', []);
xlabel('Trials from tt50');

subplot(4,4,15); hold on;
cData = perm_subtract;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc(cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(tt50_rev_perm, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);

subplot(4,4,16); hold on;
cData = aligned_subtract_perm(ix_perm, :);
imagesc('XData', xData,  'CData',  cData, clim);
plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', []);
xlabel('Trials from tt50');



formatFigurePoster([8 10], '', 10);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

% second averages
savename = ['tt50_aligned_newCsPlus_avgs' '_' expType];
ensureFigure(savename, 1);
marker = '.';
subplot(2,2,1); title('new Cs+ (acquisition)');
[hl, hp] = boundedline(xData, [nanmean(aligned_licks)' nanmean(aligned_licks_perm)'], permute([nanSEM(aligned_licks)' nanSEM(aligned_licks)'], [1 3 2]),...
    'cmap', [0 0 0; 0.7 0.7 0.7], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
set(hl, 'Marker', marker);
ylabel('cue licks');
addOrginLines;
legend(hl, 'Intact', 'Permuted', 'Location', 'best'); legend boxoff;

subplot(2,2,2);
[hl, hp] = boundedline(xData, [nanmean(aligned_ch1)' nanmean(aligned_ch1_perm)'], permute([nanSEM(aligned_ch1)' nanSEM(aligned_ch1_perm)'], [1 3 2]),...
    'cmap', [mycolors('ChAT'); 0.7 0.7 0.7], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
set(hl, 'Marker', marker);
ylabel('ACh.');
addOrginLines;

subplot(2,2,3);
[hl, hp] = boundedline(xData, [nanmean(aligned_ch2)' nanmean(aligned_ch2_perm)'], permute([nanSEM(aligned_ch2)' nanSEM(aligned_ch2_perm)'], [1 3 2]),...
    'cmap', [mycolors('DAT'); 0.7 0.7 0.7], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
set(hl, 'Marker', marker);
ylabel('Dop.');
xlabel('trials from lick changepoint');
addOrginLines;

subplot(2,2,4);
[hl, hp] = boundedline(xData, [nanmean(aligned_subtract)' nanmean(aligned_subtract_perm)'], permute([nanSEM(aligned_subtract)' nanSEM(aligned_subtract_perm)'], [1 3 2]),...
    'cmap', [0 1 0; 0.7 0.7 0.7], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
ylabel('ACh. - Dop.');
xlabel('trials from lick changepoint');
addOrginLines;

formatFigurePoster([5.5 4], '', 8);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    













%% despite paired measurements, weak correlations between detected changepoints
savename = ['ChangePoint_correlations_scatter' '_' expType];
ensureFigure(savename, 1);
subplot(3,2,1); hold on; title('new Cs+'); set(gca, 'XLim', [-20 50], 'YLim', [-20 50]); addUnityLine;
xData = all_cps(:,1); yData = all_cps(:,2);
scatter(xData, yData, 10);
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5);    
     xlabel('licks'); ylabel('ACh.');
     setXYsameLimit;
     t=textBox(sprintf('Rsq=%.3g', corr(xData, yData))); t.Interpreter='tex';t.FontSize=10;
subplot(3,2,2); hold on; title('new Cs-'); addUnityLine; % set(gca, 'XLim', [-20 50], 'YLim', [-20 50]); addUnityLine;
xData = all_cps(:,4); yData = all_cps(:,5);
scatter(xData, yData, 10); 
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5);
    xlabel('licks'); ylabel('ACh.');
    setXYsameLimit;
    t=textBox(sprintf('Rsq=%.3g', corr(xData, yData))); t.Interpreter='tex';t.FontSize=10;
subplot(3,2,3); hold on; addUnityLine; %set(gca, 'XLim', [-20 50], 'YLim', [-20 50]); addUnityLine;
xData = all_cps(:,1); yData = all_cps(:,3);
scatter(xData, yData, 10); 
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5);
    xlabel('licks'); ylabel('Dop.');
    setXYsameLimit;
    t=textBox(sprintf('Rsq=%.3g', corr(xData, yData))); t.Interpreter='tex';t.FontSize=10;
subplot(3,2,4); hold on; addUnityLine; % set(gca, 'XLim', [-20 50], 'YLim', [-20 50]); addUnityLine;
xData = all_cps(:,4); yData = all_cps(:,6);
scatter(xData, yData, 10); 
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5);
    xlabel('licks'); ylabel('Dop.');
    setXYsameLimit;
    t=textBox(sprintf('Rsq=%.3g', corr(xData, yData))); t.Interpreter='tex';t.FontSize=10;
subplot(3,2,5); hold on; addUnityLine; % set(gca, 'XLim', [-20 50], 'YLim', [-20 50]); addUnityLine;
xData = all_cps(:,2); yData = all_cps(:,3);
scatter(xData, yData, 10); 
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5);
    xlabel('Ach.'); ylabel('Dop.');
    setXYsameLimit;
    t=textBox(sprintf('Rsq=%.3g', corr(xData, yData))); t.Interpreter='tex';t.FontSize=10;
subplot(3,2,6); hold on; addUnityLine; % set(gca, 'XLim', [-20 50], 'YLim', [-20 50]); addUnityLine;
xData = all_cps(:,5); yData = all_cps(:,6);
scatter(xData, yData, 10); 
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5);
    xlabel('Ach.'); ylabel('Dop.');
    setXYsameLimit;
    t=textBox(sprintf('Rsq=%.3g', corr(xData, yData))); t.Interpreter='tex'; t.FontSize=10;
fig = gcf;
 sameXYScale(fig.Children);
fig.Position = [883   288   316   463];
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    
%     ha = findobj(gcf, 'Type', 'Axes');
%     sameXYScale(ha);



%% plot example reversals with tt50 fits, changepoints
nShow = 6;

ordering = randperm(sum(goodReversals), nShow);
gr = find(goodReversals);

toShow = gr(ordering);
% toShow = [86    70    42    26    77    44];
% nice toShow for csPlus, phPeakMean_cs_ch1: [86    70    42    26    77    44]


condition = 'csPlus'; % csPlus for newCsPlus
dataField = 'phPeakMean_cs_ch1';
fieldLabel = 'ACh.';
savename = ['reversals_pooled_Latency_examples_' fieldLabel];
ensureFigure(savename, 1);
% dataField = 'licks_cs';
for counter = 1:nShow   
    thisRev = toShow(counter);
%     thisRevCP = 
    
    % tt50
    subplot(nShow, 3, counter*3 - 2); hold on;
    if counter == 1
        title('Weibull');
    end
    plot(tt50.(condition).(dataField).toFit{thisRev}, '.', 'Color', [0 0.75 0]); 
    plot(tt50.(condition).(dataField).object{thisRev}); legend off;
    plot([1 1] + baselineTrials, get(gca, 'YLim'), '--k');
        set(gca, 'XLim', [0 120]);
    if counter == round(nShow/2)
        ylabel(fieldLabel, 'Interpreter', 'none', 'FontWeight', 'bold');
    else
        ylabel('');
    end
    % changepoint
    subplot(nShow, 3, counter*3 - 1); hold on;
    if counter == 1
        title('Changepoint');
    end    
    scatter(1:length(tt50.(condition).(dataField).cumsum{thisRev}), tt50.(condition).(dataField).cumsum{thisRev}, 10, tt50.(condition).(dataField).logitAll{thisRev}); colormap jet;
    line(repmat(tt50.(condition).(dataField).tt50(thisRev), 1, 2), get(gca, 'YLim')); 
    plot([1 1] + baselineTrials, get(gca, 'YLim'), '--k');
    set(gca, 'XLim', [0 120]);
    textBox(sprintf('Logit=%.2f', tt50.(condition).(dataField).logit(thisRev)));
    % exponential
    subplot(nShow, 3, counter*3); hold on;    
    if counter == 1
        title('Exponential');
    end    
    plot(expFit.(condition).(dataField).toFit{thisRev}, '.', 'Color', [0 0.75 0]); hold on;
    plot(expFit.(condition).(dataField).object{thisRev}); legend off;
    plot([0 0], get(gca, 'YLim'), '--k');
    set(gca, 'XLim', [0 120 - baselineTrials]);
end

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end

%% images ordered by lick changepoints
% 
% fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'licks_cs', 'csLicksROC', 'pupil_cs', 'whisk_cs'};
% titles = {'ACh', 'Dop.', 'Licks', 'Lick auROC', 'Pupil', 'Whisk'};
% clim = [-5 5];
% fh=[];
% 
% sortVariable = (cp_licks.tt50 - baselineTrials); % e.g. 20 baseline trials
% sortVariable(~goodReversals) = NaN;
% sortVariable(cp_licks.logit <= 3) = NaN;
% [sorted, sortOrder] = sort(sortVariable);
% 
% 
% saveName = 'newCsPlus_image_changepoints';
% fh(end+1) = ensureFigure(saveName, 1); colormap parula;
% cLimFactor = 3;
% xData = [min(newCsPlus.trialNumber), max(newCsPlus.trialNumber)];
% xlim = [-30 70];
% for fcounter = 1:length(fieldsToShow)
%     sfield = fieldsToShow{fcounter};
%     subplot(2,3,fcounter);
%     cData = newCsPlus.(sfield)(sortOrder, :);
%     cData = smoothdata(cData, 2, 'movmean', 5, 'omitnan');
%     imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
%     line(sorted, 1:length(sortOrder), 'Color', 'w', 'LineWidth', 2); 
%     set(gca, 'YLim', [1 find(isnan(sorted), 1) - 1]);
%     set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
%     t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
% end
% subplot(2,3,5); xlabel('Odor presentations from reversal');
% subplot(2,3,2); title('New Cs+');
% set(gcf, 'Position', [304   217   633   485]);
% 
% % %% data and images aligned to lick changepoint
% % fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'licks_cs', 'whisk_cs'};
% % titles = {'ACh', 'Dop.', 'Licks', 'Whisk'};
% % clim = [-5 5];
% % fh=[];
% % 
% % for fcounter = 1:length(fieldsToShow)
% %     sfield = fieldsToShow{fcounter};
% %     subplot(2,2,fcounter);
%     


%% cross correlations, new Cs+, try zscoring 
savename = 'xcorr_newCsPlus';
ensureFigure(savename, 1); 
% [R, lags] = avgXCorr(x, y, maxlag)
    maxlag = 20;
    allR = zeros(maxlag*2+1, sum(goodReversals), 6);
    lags = -maxlag:maxlag;
%     chat = nanzscore(newCsPlus.phPeakMean_cs_ch1, 0, 2);
%     dat = nanzscore(newCsPlus.phPeakMean_cs_ch2, 0, 2);
%     licks = nanzscore(newCsPlus.licks_cs, 0, 2);
    chat = newCsPlus.phPeakMean_cs_ch1;
    dat = newCsPlus.phPeakMean_cs_ch2;
    licks = newCsPlus.licks_cs;        
    theseOnes = find(goodReversals);
    for counter = 1:length(theseOnes)
        thisRev = theseOnes(counter);
        thisChat = chat(thisRev, isfinite(chat(thisRev, :)));
        thisDat = dat(thisRev, isfinite(dat(thisRev, :)));
        theseLicks = licks(thisRev, isfinite(licks(thisRev, :)));
        thisChatZ = zscore(thisChat);
        thisDatZ = zscore(thisDat);
        theseLicksZ = zscore(theseLicks);
        [theseR, ~] = xcorr(thisChat, thisDat, maxlag, 'coeff');
        allR(:,counter, 1) = theseR;
        [theseR, ~] = xcorr(thisChatZ, thisDatZ, maxlag, 'coeff');
        allR(:,counter, 2) = theseR;        
        [theseR, ~] = xcorr(thisChat, theseLicks, maxlag, 'coeff');
        allR(:,counter, 3) = theseR;                
        [theseR, ~] = xcorr(thisChatZ, theseLicksZ, maxlag, 'coeff');
        allR(:,counter, 4) = theseR;                
        [theseR, ~] = xcorr(thisDat, theseLicks, maxlag, 'coeff');
        allR(:,counter, 5) = theseR;                        
        [theseR, ~] = xcorr(thisDatZ, theseLicksZ, maxlag, 'coeff');
        allR(:,counter, 6) = theseR;                
    end
    R = squeeze(nanmean(allR, 2));

subplot(3,2,1); plot(lags, R(:,1)); title('chat vs dat'); 
subplot(3,2,2); plot(lags, R(:,2)); title('chat vs dat (ZScore)'); 
subplot(3,2,3); plot(lags, R(:,3)); title('chat vs licks'); 
subplot(3,2,4); plot(lags, R(:,4)); title('chat vs licks (ZScore)');
subplot(3,2,5); plot(lags, R(:,5)); title('dat vs licks'); xlabel('new cs+ trials');
subplot(3,2,6); plot(lags, R(:,6)); title('dat vs licks (ZScore)'); xlabel('new cs+ trials');


if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end


%% cross correlations, new Cs+, 
savename = 'xcorr_newCsPlus';
ensureFigure(savename, 1); 
% [R, lags] = avgXCorr(x, y, maxlag)
    maxlag = 20;
    allR = zeros(maxlag*2+1, sum(goodReversals), 6);
    lags = -maxlag:maxlag;
%     chat = nanzscore(newCsPlus.phPeakMean_cs_ch1, 0, 2);
%     dat = nanzscore(newCsPlus.phPeakMean_cs_ch2, 0, 2);
%     licks = nanzscore(newCsPlus.licks_cs, 0, 2);
    chat = newCsPlus.phPeakMean_cs_ch1;
    dat = newCsPlus.phPeakMean_cs_ch2;
    licks = newCsPlus.licks_cs;        
    theseOnes = find(goodReversals);
    for counter = 1:length(theseOnes)
        thisRev = theseOnes(counter);
        thisChat = chat(thisRev, isfinite(chat(thisRev, :)));
        thisDat = dat(thisRev, isfinite(dat(thisRev, :)));
        theseLicks = licks(thisRev, isfinite(licks(thisRev, :)));
        [theseR, ~] = xcorr(thisChat, thisDat, maxlag, 'coeff');
        allR(:,counter, 1) = theseR;
        [theseR, ~] = xcorr(thisChat, thisChat, maxlag, 'coeff');
        allR(:,counter, 2) = theseR;        
        [theseR, ~] = xcorr(thisDat, thisDat, maxlag, 'coeff');
        allR(:,counter, 3) = theseR;                
        [theseR, ~] = xcorr(thisChat, theseLicks, maxlag, 'coeff');
        allR(:,counter, 4) = theseR;                
        [theseR, ~] = xcorr(thisDat, theseLicks, maxlag, 'coeff');
        allR(:,counter, 5) = theseR;                        
        [theseR, ~] = xcorr(theseLicks, theseLicks, maxlag, 'coeff');
        allR(:,counter, 6) = theseR;                
    end
    R = squeeze(nanmean(allR, 2));

subplot(3,2,1); plot(lags, R(:,1)); title('chat vs dat'); 
subplot(3,2,2); plot(lags, R(:,2)); title('chat'); 
subplot(3,2,3); plot(lags, R(:,3)); title('dat'); 
subplot(3,2,4); plot(lags, R(:,4)); title('chat vs licks');
subplot(3,2,5); plot(lags, R(:,5)); title('dat vs licks'); xlabel('new cs+ trials');
subplot(3,2,6); plot(lags, R(:,6)); title('licks'); xlabel('new cs+ trials');
fig=gcf;
set(fig.Children, 'YLim', [0 1]);


if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end


%% cross correlations, new Cs+, corrected

savename = 'xcorr_newCsPlus_corrected';
ensureFigure(savename, 1); 
% [R, lags] = avgXCorr(x, y, maxlag)
maxlag = 10;
trialRange = [30 30];
R = deal(zeros(maxlag*2+1, 6, 3)); % corrected, shift predictor, raw occupy third dimension
lags = -maxlag:maxlag;

%     chat = nanzscore(newCsPlus.phPeakMean_cs_ch1, 0, 2);
%     dat = nanzscore(newCsPlus.phPeakMean_cs_ch2, 0, 2);
%     licks = nanzscore(newCsPlus.licks_cs, 0, 2);

chat = newCsPlus.phPeakMean_cs_ch1(goodReversals, newCsPlus.firstRevTrial - trialRange:newCsPlus.firstRevTrial + trialRange - 1);
dat = newCsPlus.phPeakMean_cs_ch2(goodReversals, newCsPlus.firstRevTrial - trialRange:newCsPlus.firstRevTrial + trialRange - 1);
licks = newCsPlus.licks_cs(goodReversals, newCsPlus.firstRevTrial - trialRange:newCsPlus.firstRevTrial + trialRange - 1);        

[R1, R2, R3, testlags] = correctedXCorr(chat, dat, maxlag, 2);
R(:,1,1) = R1; R(:,1,2) = R2; R(:,1,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(chat, chat, maxlag, 2);
R(:,2,1) = R1; R(:,2,2) = R2; R(:,2,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(dat, dat, maxlag, 2);
R(:,3,1) = R1; R(:,3,2) = R2; R(:,3,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(chat, licks, maxlag, 2);
R(:,4,1) = R1; R(:,4,2) = R2; R(:,4,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(dat, licks, maxlag, 2);
R(:,5,1) = R1; R(:,5,2) = R2; R(:,5,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(licks, licks, maxlag, 2);
R(:,6,1) = R1; R(:,6,2) = R2; R(:,6,3) = R3;
    
   

subplot(3,2,1); plot(lags, squeeze(R(:,1,:))); title('chat vs dat'); legend('corrected', 'shift predictor', 'raw'); legend boxoff;
subplot(3,2,2); plot(lags, squeeze(R(:,2,:))); title('chat');  legend('corrected', 'shift predictor', 'raw'); legend boxoff;
subplot(3,2,3); plot(lags, squeeze(R(:,3,:))); title('dat'); 
subplot(3,2,4); plot(lags, squeeze(R(:,4,:))); title('chat vs licks');
subplot(3,2,5); plot(lags, squeeze(R(:,5,:))); title('dat vs licks'); xlabel('new cs+ trials');
subplot(3,2,6); plot(lags, squeeze(R(:,6,:))); title('licks'); xlabel('new cs+ trials');


if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end


%% cross correlations, new Cs-, corrected

savename = 'xcorr_newCsMinus_corrected';
ensureFigure(savename, 1); 
% [R, lags] = avgXCorr(x, y, maxlag)
maxlag = 10;
trialRange = [30 30];
R = deal(zeros(maxlag*2+1, 6, 3)); % corrected, shift predictor, raw occupy third dimension
lags = -maxlag:maxlag;

%     chat = nanzscore(newCsPlus.phPeakMean_cs_ch1, 0, 2);
%     dat = nanzscore(newCsPlus.phPeakMean_cs_ch2, 0, 2);
%     licks = nanzscore(newCsPlus.licks_cs, 0, 2);

chat = newCsMinus.phPeakMean_cs_ch1(goodReversals, newCsMinus.firstRevTrial - trialRange:newCsMinus.firstRevTrial + trialRange - 1);
dat = newCsMinus.phPeakMean_cs_ch2(goodReversals, newCsMinus.firstRevTrial - trialRange:newCsMinus.firstRevTrial + trialRange - 1);
licks = newCsMinus.licks_cs(goodReversals, newCsMinus.firstRevTrial - trialRange:newCsMinus.firstRevTrial + trialRange - 1);        

[R1, R2, R3, testlags] = correctedXCorr(chat, dat, maxlag, 2);
R(:,1,1) = R1; R(:,1,2) = R2; R(:,1,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(chat, chat, maxlag, 2);
R(:,2,1) = R1; R(:,2,2) = R2; R(:,2,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(dat, dat, maxlag, 2);
R(:,3,1) = R1; R(:,3,2) = R2; R(:,3,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(chat, licks, maxlag, 2);
R(:,4,1) = R1; R(:,4,2) = R2; R(:,4,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(dat, licks, maxlag, 2);
R(:,5,1) = R1; R(:,5,2) = R2; R(:,5,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(licks, licks, maxlag, 2);
R(:,6,1) = R1; R(:,6,2) = R2; R(:,6,3) = R3;
    
   

subplot(3,2,1); plot(lags, squeeze(R(:,1,:))); title('chat vs dat'); legend('corrected', 'shift predictor', 'raw'); legend boxoff;
subplot(3,2,2); plot(lags, squeeze(R(:,2,:))); title('chat');  legend('corrected', 'shift predictor', 'raw'); legend boxoff;
subplot(3,2,3); plot(lags, squeeze(R(:,3,:))); title('dat'); 
subplot(3,2,4); plot(lags, squeeze(R(:,4,:))); title('chat vs licks');
subplot(3,2,5); plot(lags, squeeze(R(:,5,:))); title('dat vs licks'); xlabel('new cs- trials');
subplot(3,2,6); plot(lags, squeeze(R(:,6,:))); title('licks'); xlabel('new cs- trials');


if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end

%% cross correlations, new Cs-, 
savename = 'xcorr_newCsMinus';
ensureFigure(savename, 1); 
% [R, lags] = avgXCorr(x, y, maxlag)
    maxlag = 20;
    allR = zeros(maxlag*2+1, sum(goodReversals), 6);
    lags = -maxlag:maxlag;
%     chat = nanzscore(newCsMinus.phPeakMean_cs_ch1, 0, 2);
%     dat = nanzscore(newCsMinus.phPeakMean_cs_ch2, 0, 2);
%     licks = nanzscore(newCsMinus.licks_cs, 0, 2);
    chat = newCsMinus.phPeakMean_cs_ch1;
    dat = newCsMinus.phPeakMean_cs_ch2;
    licks = newCsMinus.licks_cs;    
    theseOnes = find(goodReversals);
    for counter = 1:length(theseOnes)
        thisRev = theseOnes(counter);
        thisChat = chat(thisRev, isfinite(chat(thisRev, :)));
        thisDat = dat(thisRev, isfinite(dat(thisRev, :)));
        theseLicks = licks(thisRev, isfinite(licks(thisRev, :)));
        [theseR, ~] = xcorr(thisChat, thisDat, maxlag, 'coeff');
        allR(:,counter, 1) = theseR;
        [theseR, ~] = xcorr(thisChat, thisChat, maxlag, 'coeff');
        allR(:,counter, 2) = theseR;        
        [theseR, ~] = xcorr(thisDat, thisDat, maxlag, 'coeff');
        allR(:,counter, 3) = theseR;                
        [theseR, ~] = xcorr(thisChat, theseLicks, maxlag, 'coeff');
        allR(:,counter, 4) = theseR;                
        [theseR, ~] = xcorr(thisDat, theseLicks, maxlag, 'coeff');
        allR(:,counter, 5) = theseR;                        
        [theseR, ~] = xcorr(theseLicks, theseLicks, maxlag, 'coeff');
        allR(:,counter, 6) = theseR;                
    end
    R = squeeze(nanmean(allR, 2));

subplot(3,2,1); plot(lags, R(:,1)); title('chat vs dat'); 
subplot(3,2,2); plot(lags, R(:,2)); title('chat'); 
subplot(3,2,3); plot(lags, R(:,3)); title('dat'); 
subplot(3,2,4); plot(lags, R(:,4)); title('chat vs licks');
subplot(3,2,5); plot(lags, R(:,5)); title('dat vs licks'); xlabel('new cs- trials');
subplot(3,2,6); plot(lags, R(:,6)); title('licks'); xlabel('new cs- trials');


if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end

%% test correlations
x = linspace(0, 1, 10000);
y1 = 1 * x + 40;
y2 =  1 * x;

savename = 'XCorr_baseline_effect';
ensureFigure(savename, 1);
subplot(2,2,1); plot(x', [y1' y2']); title('not zscored'); xlabel('x'); ylabel('y'); set(gca, 'YLim', [-10 50]);
[r,lags] = xcorr(y1,y2, 20, 'unbiased');
subplot(2,2,3); plot(lags, r); ylabel('xcorr'); xlabel('lags');
subplot(2,2,2); plot(x ,zscore(y1), '-', 'LineWidth', 2); hold on; plot(x ,zscore(y2), '--', 'LineWidth', 4);  title('zscored'); xlabel('x'); ylabel('y (zscored)');
[r,lags] = xcorr(y1 - mean(y1) ,y2 - mean(y2), 20, 'unbiased');
subplot(2,2,4); plot(lags, r); ylabel('xcorr'); xlabel('lags');

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end
%% compare strength of licking correlations, ChAT-cre vs DAT-cre
mouseNumber_newCsPlus = repmat(mouseNumber, 1, size(newCsPlus.licks_cs, 2));
mouseNumber_newCsMinus = repmat(mouseNumber, 1, size(newCsMinus.licks_cs, 2));
all_Licks = [reshape(newCsPlus.licks_cs(goodReversals, :), numel(newCsPlus.licks_cs(goodReversals, :)), 1); reshape(newCsMinus.licks_cs(goodReversals, :), numel(newCsMinus.licks_cs(goodReversals, :)), 1)]; 
all_ChAT = [reshape(newCsPlus.phPeakMean_cs_ch1(goodReversals, :), numel(newCsPlus.phPeakMean_cs_ch1(goodReversals, :)), 1); reshape(newCsMinus.phPeakMean_cs_ch1(goodReversals, :), numel(newCsMinus.phPeakMean_cs_ch1(goodReversals, :)), 1)]; 
all_DAT = [reshape(newCsPlus.phPeakMean_cs_ch2(goodReversals, :), numel(newCsPlus.phPeakMean_cs_ch2(goodReversals, :)), 1); reshape(newCsMinus.phPeakMean_cs_ch2(goodReversals, :), numel(newCsMinus.phPeakMean_cs_ch2(goodReversals, :)), 1)]; 
all_ChAT_us = [reshape(newCsPlus.phPeakMean_us_ch1(goodReversals, :), numel(newCsPlus.phPeakMean_us_ch1(goodReversals, :)), 1); reshape(newCsMinus.phPeakMean_us_ch1(goodReversals, :), numel(newCsMinus.phPeakMean_us_ch1(goodReversals, :)), 1)]; 
all_DAT_us = [reshape(newCsPlus.phPeakMean_us_ch2(goodReversals, :), numel(newCsPlus.phPeakMean_us_ch2(goodReversals, :)), 1); reshape(newCsMinus.phPeakMean_us_ch2(goodReversals, :), numel(newCsMinus.phPeakMean_us_ch2(goodReversals, :)), 1)]; 
all_mouseNumber = [reshape(mouseNumber_newCsPlus(goodReversals, :), numel(mouseNumber_newCsPlus(goodReversals, :)), 1); reshape(mouseNumber_newCsMinus(goodReversals, :), numel(mouseNumber_newCsMinus(goodReversals, :)), 1)]; 
keep = isfinite(all_Licks) & isfinite(all_ChAT) & isfinite(all_DAT);
all_Licks = all_Licks(keep);
all_ChAT = all_ChAT(keep);
all_DAT = all_DAT(keep);

all_ChAT_us = all_ChAT_us(keep);
all_DAT_us = all_DAT_us(keep);
all_mouseNumber = all_mouseNumber(keep);

ensureFigure('test_corr', 1); 
subplot(1,2,1); scatter(all_Licks, all_ChAT, 8, '.'); ylabel('cue ChAT'); textBox(sprintf('R=%.2f', corr(all_Licks, all_ChAT)));
subplot(1,2,2); scatter(all_Licks, all_DAT, 8, '.'); ylabel('cue DAT'); textBox(sprintf('R=%.2f', corr(all_Licks, all_DAT)));

%% cs and us correlations between trials


linecolors = [1 0 0; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 1];     
keepers = [1:6];

nMice = max(mouseNumber);
savename = 'reversal_trial_correlations';
ensureFigure(savename, 1);
% allTrials = true; 
ll = [];
subplot(1,2,1); hold on; subplot(1,2,2); hold on;
for counter = 1:nMice
%     if ~ismember(counter, keepers)
%         continue;
%     end
    thisMouse = all_mouseNumber == counter;
    subplot(1,2,1); %set(gca, 'XLim', xlims(1,:));
%     allTrials = allTrials | trialSets{counter};
    xData = all_ChAT(thisMouse); yData = all_DAT(thisMouse); 
%     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
    scatter(xData, yData, 12, linecolors(counter, :), '.');
    % fit for cs
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob, 'predfunc'); legend off; %,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
    ll(end + 1) = fph(1);
    
    subplot(1,2,2); %set(gca, 'XLim', xlims(1,:));
%     allTrials = allTrials | trialSets{counter};
    xData = all_ChAT_us(thisMouse); yData = all_DAT_us(thisMouse); 
%     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for cs
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob, 'predfunc'); legend off; %,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));    
end

subplot(1,2,1); legend(ll, DB.animals(keepers), 'Interpreter', 'none');
title('Cs'); ylabel('Dop.'); xlabel('ACh.');
subplot(1,2,2);
title('Us'); ylabel('Dop.'); xlabel('ACh.');

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end

%% example of changepoint detection


savename = 'changepoint_explanation';
ensureFigure(savename, 1);

thisGoodRev = 4;
goodOnes2 = find(goodOnes);
thisRev = goodOnes2(thisGoodRev);
subplot(1,2,1); hold on;

rawData = [AR.csMinus.licks_cs.before(thisRev, end - baselineTrials + 1:end) AR.csPlus.licks_cs.after(thisRev, 1:end)];
valid = find(~isnan(rawData));
rawData = rawData(valid);
plot(rawData, '-b');
plot([baselineTrials baselineTrials] + 1, get(gca, 'YLim'), '--k');
title('Cue lick data');
set(gca, 'XLim', [0 length(valid)]);
ylabel('licks');
xlabel('trials spanning reversal');
subplot(1,2,2); hold on;
title('Changepoint detection');



np = length(rawData);
total = sum(rawData);
nochange = total/np * (1:np); nochange = nochange(:); % diagonal line from orgin

plot(nochange, '-', 'Color', [0.7 0.7 0.7]); ylabel('cumulative licks');
scatter(1:length(tt50.csPlus.licks_cs.cumsum{thisRev}), tt50.csPlus.licks_cs.cumsum{thisRev}, 10, tt50.csPlus.licks_cs.logitAll{thisRev}); colormap jet;

line(repmat(tt50.csPlus.licks_cs.tt50(thisRev), 1, 2), get(gca, 'YLim'), 'Color', [0 1 0]); 
plot([1 1] + baselineTrials, get(gca, 'YLim'), '--k');
set(gca, 'XLim', [0 120]);
t = textBox(sprintf('Logit=%.2f', tt50.csPlus.licks_cs.logit(thisRev)));
t.Color = [0 0.5 0];

set(gca, 'XLim', [0 length(valid)]);
xlabel('trials spanning reversal');

formatFigurePoster([6 3], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end   


%%
% 
% minLogit = 2;
% zeroTrials = cp.csPlus.licks_cs.tt50;
% zeroTrials(cp.csPlus.licks_cs.logit < minLogit) = NaN;
% zeroTrials = zeroTrials - baselineTrials;
% trialWindow = [-20 50];
% reversalPoints = 0 - zeroTrials;
% [sorted, ix] = sort(reversalPoints);
% [aligned_subtract, xData] = alignedDataWindow(newCsPlus.phPeakMean_cs_AchMinusDop, true(nReversals,1), 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1), nReversals, 1));
% [aligned_licks, ~] = alignedDataWindow(newCsPlus.licks_cs, true(nReversals,1), 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1), nReversals, 1));
% [aligned_ch1, ~] = alignedDataWindow(newCsPlus.phPeakMean_cs_ch1, true(nReversals,1), 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1), nReversals, 1));
% [aligned_ch2, ~] = alignedDataWindow(newCsPlus.phPeakMean_cs_ch2, true(nReversals,1), 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1), nReversals, 1));
% 
% 
% savename = 'cp_aligned_newCsPlus_images';
% ensureFigure(savename, 1);
% 
% subplot(2,2,1); hold on;
% % first images, sorted by reversal point.
% title('Cue licks');
% imagesc('XData', xData, 'CData', aligned_licks(goodReversals & isfinite(zeroTrials), :));
% plot(reversalPoints(ix(goodReversals & isfinite(zeroTrials))), 1:sum(goodReversals & isfinite(zeroTrials)), '--w', 'LineWidth', 2); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodReversals & isfinite(zeroTrials))]);
% 
% subplot(2,2,2); hold on;
% title('ACh.', 'Color', mycolors('ChAT'));
% imagesc('XData', xData, 'CData', aligned_ch1(ix(goodReversals & isfinite(zeroTrials)), :));
% plot(reversalPoints(ix(goodReversals & isfinite(zeroTrials))), 1:sum(goodReversals & isfinite(zeroTrials)), '--w', 'LineWidth', 2); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodReversals & isfinite(zeroTrials))]);
% 
% subplot(2,2,3); hold on;
% title('Dop.', 'Color', mycolors('DAT'));
% imagesc('XData', xData, 'CData', aligned_ch2(ix(goodReversals & isfinite(zeroTrials)), :));
% plot(reversalPoints(ix(goodReversals & isfinite(zeroTrials))), 1:sum(goodReversals & isfinite(zeroTrials)), '--w', 'LineWidth', 2); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodReversals & isfinite(zeroTrials))]);
% 
% subplot(2,2,4); hold on;
% title('ACh. - Dop.');
% imagesc('XData', xData, 'CData', aligned_subtract(ix(goodReversals & isfinite(zeroTrials)), :));
% plot(reversalPoints(ix(g

