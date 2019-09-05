DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);

saveOn = 1;


%%
exp.value = {'DC_44'  'DC_46'  'DC_47' 'DC_53'  'DC_54'  'DC_56'}; % exclude DC_51
exp.valence = {'DC_17'  'DC_20'  'DC_35'  'DC_36'  'DC_37'  'DC_40'};
exp.all = [exp.value exp.valence];

%% averages
expType = 'valence';
smoothWindow = 3;
compile_reversal_data;


trialWindow = [-30 30];
ylim = [-1 2];
figsize = [1.4 1.2];
% new cs plus
% common = sum(~isnan(newCsPlus.licks_cs)) > 3;
common = newCsPlus.firstRevTrial + trialWindow(1):newCsPlus.firstRevTrial + trialWindow(2) - 1;
savename = 'reversals_newCsPlus';
ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsPlus.trialNumber(common), nanmean(newCsPlus.licks_cs(goodReversals, common)), nanSEM(newCsPlus.licks_cs(goodReversals, common))',...
    'cmap', mycolors('licks'), 'nan', 'gap');
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus.trialNumber(common), nanmean(newCsPlus.phPeakMean_cs_ch2(goodReversals, common)), nanSEM(newCsPlus.phPeakMean_cs_ch2(goodReversals, common))',...
    'cmap', mycolors('dat'), 'nan', 'gap'); hold on
hla(2) = hl;
[hl, hp] = boundedline(newCsPlus.trialNumber(common), nanmean(newCsPlus.phPeakMean_cs_ch1(goodReversals, common)), nanSEM(newCsPlus.phPeakMean_cs_ch1(goodReversals, common))',...
    'cmap', mycolors('chat'), 'nan', 'gap');
hla(3) = hl;


set(hla, 'LineWidth', 1);
set(gca, 'XLim', trialWindow);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 1);
legend(hla, { '\bf\color[rgb]{0.5,0.5,0.5}Licks', '\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.'},...
    'Location', 'northwest', 'FontSize', 6, 'Interpreter', 'tex', 'Box', 'off');

title('New Cs+');
xlabel('Odor presentations from rev.');
ylabel('Cue response');

formatFigurePublish('size', figsize, 'fontSize', 6);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
end   


% new cs minus
% common = sum(~isnan(newCsMinus.licks_cs)) > 3;
common = newCsMinus.firstRevTrial + trialWindow(1):newCsMinus.firstRevTrial + trialWindow(2) - 1;

savename = 'reversals_newCsMinus';
ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsMinus.trialNumber(common), nanmean(newCsMinus.licks_cs(goodReversals, common)), nanSEM(newCsMinus.licks_cs(goodReversals, common))',...
    'cmap', mycolors('licks'));
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus.trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_ch2(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_ch2(goodReversals, common))',...
    'cmap', mycolors('dat')); hold on
hla(2) = hl;
[hl, hp] = boundedline(newCsMinus.trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_ch1(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_ch1(goodReversals, common))',...
    'cmap', mycolors('chat'));
hla(3) = hl;

set(hla, 'LineWidth', 1);
set(gca, 'XLim', trialWindow);%, 'YLim', ylim);    
h  = addOrginLines;
set(h, 'LineWidth', 1);
legend(hla, {'\bf\color[rgb]{0.5,0.5,0.5}Licks', '\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.'},...
            'Location', 'northeast', 'FontSize', 6, 'Interpreter', 'tex', 'Box', 'off');
title('New Cs-');
xlabel('Odor presentations from rev.');
ylabel('Cue response');    
formatFigurePublish('size', figsize, 'fontSize', 6);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
end   


%% subtraction averages- show run-down and difference upon extinction

trialWindow = [-30 30];
ylim = [-1 2];
figsize = [1.4 1.2];
% new cs plus
% common = sum(~isnan(newCsPlus.licks_cs)) > 3;


savename = ['ACh_minus_Dop_avgs_newCsPlus' '_' expType];
ensureFigure(savename, 1);
common = newCsPlus.firstRevTrial + trialWindow(1):newCsPlus.firstRevTrial + trialWindow(2) - 1;
axes;
title('new Cs+');
[hl, hp] = boundedline(newCsPlus.trialNumber(common), nanmean(newCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common)), nanSEM(newCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common))',...
    'cmap', [0 1 0], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
addOrginLines;
ylabel('ACh. - Dop.');    
xlabel('Odor presentations from rev.');
formatFigurePublish('size', figsize, 'fontSize', 6);

if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
end   


savename = ['ACh_minus_Dop_avgs_newCsMinus' '_' expType];
ensureFigure(savename, 1);
common = newCsMinus.firstRevTrial + trialWindow(1):newCsMinus.firstRevTrial + trialWindow(2) - 1;
axes;
title('new Cs-');
[hl, hp] = boundedline(newCsMinus.trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_AchMinusDop(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_AchMinusDop(goodReversals, common))',...
    'cmap', [0 1 0], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
addOrginLines;
ylabel('ACh. - Dop.');  
xlabel('Odor presentations from rev.');
formatFigurePublish('size', figsize, 'fontSize', 6);

if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
end   

%% changepoint detection for acquisition (acq) and extinction (ext)

expType = 'all';
smoothWindow = 1; % smoothing ruins cross correlation validity, also problematic for permutation test
compile_reversal_data;
minLogit = 3;

baselineTrials = 20;
cp.csPlus.licks_cs = bpChangePoints(newCsPlus.licks_cs(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
cp.csPlus.phPeakMean_cs_ch1 = bpChangePoints(newCsPlus.phPeakMean_cs_ch1(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
cp.csPlus.phPeakMean_cs_ch2 = bpChangePoints(newCsPlus.phPeakMean_cs_ch2(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
cp.csMinus.licks_cs = bpChangePoints(newCsMinus.licks_cs(:,newCsMinus.firstRevTrial - baselineTrials:end), 2, 1000, 'down');
cp.csMinus.phPeakMean_cs_ch1 = bpChangePoints(newCsMinus.phPeakMean_cs_ch1(:,newCsMinus.firstRevTrial - baselineTrials:end), 2, 1000, 'down');
cp.csMinus.phPeakMean_cs_ch2 = bpChangePoints(newCsMinus.phPeakMean_cs_ch2(:,newCsMinus.firstRevTrial - baselineTrials:end), 2, 1000, 'down');

goodPlus = goodReversals & (cp.csPlus.licks_cs.logit > minLogit) & (cp.csPlus.phPeakMean_cs_ch1.logit > minLogit) & (cp.csPlus.phPeakMean_cs_ch2.logit > minLogit);
goodMinus = goodReversals & (cp.csMinus.licks_cs.logit > minLogit) & (cp.csMinus.phPeakMean_cs_ch1.logit > minLogit) & (cp.csMinus.phPeakMean_cs_ch2.logit > minLogit);

%% make a bar graph or box plot of changepoints
figSize = [1.6 1.2];    
savename = ['cp_violin' '_' expType];
ensureFigure(savename, 1); axes; hold on;
markerSize = 4;


fields = {...
    'csPlus', 'licks_cs', mycolors('licks'), 'licks';...
    'csPlus', 'phPeakMean_cs_ch1', mycolors('chat'), 'ACh.';...
    'csPlus', 'phPeakMean_cs_ch2', mycolors('dat'), 'Dop.';...
    'csMinus', 'licks_cs', mycolors('licks'), 'licks';...
    'csMinus', 'phPeakMean_cs_ch1', mycolors('chat'), 'ACh.';...
    'csMinus', 'phPeakMean_cs_ch2', mycolors('dat'), 'Dop.';...    
    };


all_cps = cell(size(fields, 1), 1);
violins = [];
for counter = 1:size(fields, 1)
    switch fields{counter, 1}
        case 'csPlus'
            yData = cp.(fields{counter, 1}).(fields{counter, 2}).index(goodPlus) - baselineTrials;       
        case 'csMinus'
            yData = cp.(fields{counter, 1}).(fields{counter, 2}).index(goodMinus) - baselineTrials;       
    end
    xData = zeros(size(yData)) + counter;
    all_cps{counter} = yData;    
    violin = Violin(yData, counter, 'ShowMean', false, 'ShowNotches', false, 'EdgeColor', [1 1 1], 'ViolinColor', fields{counter, 3});
    % customize appearance
    violin.ScatterPlot.delete;
    violin.ScatterPlot = scatter(xData, yData, markerSize, fields{counter, 3}, 'filled');    
    violin.MedianPlot.delete;
    violin.WhiskerPlot.delete
    violin.BoxPlot.delete
    violin.ViolinAlpha = 1;
    line([counter - 0.4 counter - 0.4 counter + 0.4 counter + 0.4], [0 median(yData) median(yData) 0], 'Color', fields{counter, 3});
end

ylabel('changepoint (trials from rev.)');
% bh.BaseLine.LineStyle = '--';
% bh.BaseLine.Color = [0.7 0.7 0.7];
% bh.BaseLine.LineWidth = 2;

set(gca, 'XLim', [0.5 6.5]);
plot(get(gca, 'XLim'), [0 0], '--', 'Color', [0.7 0.7 0.7]);
set(gca, 'YLim', [-baselineTrials 50], 'XTick', [1 2 3 4 5 6], 'XTickLabel', fields(:,4));
formatFigurePublish('size', figSize, 'fontSize', 6);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
%     export_fig(fullfile(savepath, savename), '-eps');
end

% stats on changepoints


comp = {'Lick_vs_ACh'; 'Lick_vs_Dop'; 'ACh_vs_Dop'};
acq_p = zeros(3,1);
acq_n = zeros(3,1);
ext_p = zeros(3,1);
ext_n = zeros(3,1);
cp_stats = table(comp, acq_p, acq_n, ext_p, ext_n);
cp_stats.acq_p(1) = signrank(cp.csPlus.licks_cs.index(goodPlus), cp.csPlus.phPeakMean_cs_ch1.index(goodPlus));
cp_stats.acq_p(2) = signrank(cp.csPlus.licks_cs.index(goodPlus), cp.csPlus.phPeakMean_cs_ch2.index(goodPlus));
cp_stats.acq_p(3) = signrank(cp.csPlus.phPeakMean_cs_ch1.index(goodPlus), cp.csPlus.phPeakMean_cs_ch2.index(goodPlus));
cp_stats.acq_n(:) = sum(goodPlus);

cp_stats.ext_p(1) = signrank(cp.csMinus.licks_cs.index(goodMinus), cp.csMinus.phPeakMean_cs_ch1.index(goodMinus));
cp_stats.ext_p(2) = signrank(cp.csMinus.licks_cs.index(goodMinus), cp.csMinus.phPeakMean_cs_ch2.index(goodMinus));
cp_stats.ext_p(3) = signrank(cp.csMinus.phPeakMean_cs_ch1.index(goodMinus), cp.csMinus.phPeakMean_cs_ch2.index(goodMinus));
cp_stats.ext_n(:) = sum(goodMinus);


if saveOn
    save(fullfile(savepath, ['cp_stats_' expType]), 'cp_stats');    
end

cp_stats




%% align reversals by lick changepoint

% csPlus
minLogit = 3;
trialWindow = [-30 30];
cpField = 'licks_cs';
cLimFactor = 3;
trialWindow_rev = [-10 30];
markerSize = 3;
lineWidth = 2;
avgfigSize = [1.2 1.2];
figSize = [1.6 1.2];


% first new csPlus (acquisition)

% setup for aligned by changepoint
goodOnes = goodReversals & (cp.csPlus.(cpField).logit > minLogit) & (cp.csPlus.(cpField).index > (baselineTrials - 0));
zeroTrials = cp.csPlus.(cpField).index(goodOnes) - baselineTrials;
reversalPoints = 0 - zeroTrials;
reversalPoints = reversalPoints;
[sorted, ix] = sort(reversalPoints); % THEN sort them


% first select good subset
good_licks = newCsPlus.licks_cs(goodOnes, :);
% % try normalizing licks
% good_licks = good_licks ./ percentile(good_licks(:,newCsPlus.firstRevTrial:end), 0.9, 2);
good_licks = nanzscore(good_licks, 0, 2);
good_ch1 = nanzscore(newCsPlus.phPeakMean_cs_ch1(goodOnes, :), 0, 2);
good_ch2 = nanzscore(newCsPlus.phPeakMean_cs_ch2(goodOnes, :), 0, 2);

reversals = true(sum(goodOnes), 1);
[aligned_licks, xData] = alignedDataWindow(good_licks, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch1, ~] = alignedDataWindow(good_ch1, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch2, ~] = alignedDataWindow(good_ch2, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));

% % randomly permute across reversals (not time) to preserve average but
% % scramble potential changepoints
% maxTrials = length(newCsPlus.trialNumber);
% 
% % generate subscript indices for permutation
% row_ix = repmat(1:maxTrials, sum(goodOnes), 1);
% col_ix = zeros(sum(goodOnes), maxTrials);
% for counter = 1:maxTrials
%     col_ix(:,counter) = randperm(sum(goodOnes))';
% end
% 
% lin_ix = sub2ind([sum(goodOnes) maxTrials], col_ix, row_ix);
% perm_licks = good_licks(lin_ix);
% perm_ch1 = good_ch1(lin_ix);
% perm_ch2 = good_ch2(lin_ix);

% cp_perm.csPlus.licks_cs = bpChangePoints(perm_licks(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
% zeroTrials_perm = cp_perm.csPlus.licks_cs.index - baselineTrials;
% reversalPoints_perm = 0 - zeroTrials_perm;
% reversalPoints_perm = reversalPoints_perm;
% [sorted_perm, ix_perm] = sort(reversalPoints_perm); % THEN sort them

% [aligned_licks_perm, ~] = alignedDataWindow(perm_licks, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
% [aligned_ch1_perm, ~] = alignedDataWindow(perm_ch1, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
% [aligned_ch2_perm, ~] = alignedDataWindow(perm_ch2, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));


% setup for alignment by reversal
cp_rev = zeroTrials;
% cp_rev_perm = zeroTrials_perm;

% first images, sorted by reversal point.
savename = ['cp_aligned_images_newCsPlus_' expType];
ensureFigure(savename, 1);

params = struct();    
params.matpos = [0 0 1 1];    
params.figmargin = [0.15 0 0.1 0.2];
params.cellmargin = [0.025 0.025 0.025 0.025];    

hax = axesmatrix(3,2,1:6,params);

axes(hax(1)); hold on;
ylabel('Cue licks');
cData = good_licks;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), markerSize, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

axes(hax(2)); hold on;
cData = aligned_licks(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':r', 'LineWidth', lineWidth);
set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

axes(hax(3)); hold on;
ylabel('ACh.', 'Color', mycolors('ChAT'));
cData = good_ch1;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), markerSize, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

axes(hax(4)); hold on;
cData = aligned_ch1(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':r', 'LineWidth', lineWidth); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

axes(hax(5)); hold on;
ylabel('Dop.', 'Color', mycolors('DAT'));
cData = good_ch2;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), markerSize, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', [0 10 20]);
xlabel('Trials from rev.');

axes(hax(6)); hold on;
cData = aligned_ch2(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':r', 'LineWidth', lineWidth); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', [-10 0 10]);
xlabel('Trials from c.p.');

formatFigurePublish('size', figSize, 'fontSize', 6);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
    print(gcf, '-dpdf', fullfile(savepath, [savename '.pdf']));
%     export_fig(fullfile(savepath, savename), '-eps');
end

%  averages
savename = ['cp_aligned_newCsPlus_avgs' '_' expType];
ensureFigure(savename, 1);
marker = '.';
yLim = [-0.8 0.8];

axes; hold on;


[hl, hp] = boundedline(xData, [nanmean(aligned_ch1)' nanmean(aligned_ch2)'], permute([nanSEM(aligned_ch1)' nanSEM(aligned_ch2)'], [1 3 2]),...
    'cmap', [mycolors('chat'); mycolors('dat')], 'nan', 'gap');
set(gca, 'YLim', yLim);
ylabel('Cue response (z-score)');
xlabel('Trials from lick c.p.');

formatFigurePublish('size', avgfigSize, 'fontSize', 6);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
%     export_fig(fullfile(savepath, savename), '-eps');
end

% now csMinus

% csMinus
minLogit = 3;
trialWindow = [-20 20];
cpField = 'licks_cs';
cLimFactor = 3;
trialWindow_rev = [-10 30];
markerSize = 3;
lineWidth = 2;



% first new csMinus (acquisition)

% setup for aligned by changepoint
goodOnes = goodReversals & (cp.csMinus.(cpField).logit > minLogit) & (cp.csMinus.(cpField).index > (baselineTrials - 0));
zeroTrials = cp.csMinus.(cpField).index(goodOnes) - baselineTrials;
reversalPoints = 0 - zeroTrials;
reversalPoints = reversalPoints;
[sorted, ix] = sort(reversalPoints); % THEN sort them


% first select good subset
good_licks = newCsMinus.licks_cs(goodOnes, :);
% % try normalizing licks
% good_licks = good_licks ./ percentile(good_licks(:,newCsMinus.firstRevTrial:end), 0.9, 2);
good_licks = nanzscore(good_licks, 0, 2);
good_ch1 = nanzscore(newCsMinus.phPeakMean_cs_ch1(goodOnes, :), 0, 2);
good_ch2 = nanzscore(newCsMinus.phPeakMean_cs_ch2(goodOnes, :), 0, 2);

reversals = true(sum(goodOnes), 1);
[aligned_licks, xData] = alignedDataWindow(good_licks, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsMinus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch1, ~] = alignedDataWindow(good_ch1, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsMinus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch2, ~] = alignedDataWindow(good_ch2, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsMinus.trialNumber(1),  sum(goodOnes), 1));

% % randomly permute across reversals (not time) to preserve average but
% % scramble potential changepoints
% maxTrials = length(newCsMinus.trialNumber);
% 
% % generate subscript indices for permutation
% row_ix = repmat(1:maxTrials, sum(goodOnes), 1);
% col_ix = zeros(sum(goodOnes), maxTrials);
% for counter = 1:maxTrials
%     col_ix(:,counter) = randperm(sum(goodOnes))';
% end
% 
% lin_ix = sub2ind([sum(goodOnes) maxTrials], col_ix, row_ix);
% perm_licks = good_licks(lin_ix);
% perm_ch1 = good_ch1(lin_ix);
% perm_ch2 = good_ch2(lin_ix);

% cp_perm.csMinus.licks_cs = bpChangePoints(perm_licks(:,newCsMinus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
% zeroTrials_perm = cp_perm.csMinus.licks_cs.index - baselineTrials;
% reversalPoints_perm = 0 - zeroTrials_perm;
% reversalPoints_perm = reversalPoints_perm;
% [sorted_perm, ix_perm] = sort(reversalPoints_perm); % THEN sort them

% [aligned_licks_perm, ~] = alignedDataWindow(perm_licks, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsMinus.trialNumber(1),  sum(goodOnes), 1));
% [aligned_ch1_perm, ~] = alignedDataWindow(perm_ch1, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsMinus.trialNumber(1),  sum(goodOnes), 1));
% [aligned_ch2_perm, ~] = alignedDataWindow(perm_ch2, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsMinus.trialNumber(1),  sum(goodOnes), 1));


% setup for alignment by reversal
cp_rev = zeroTrials;
% cp_rev_perm = zeroTrials_perm;

% first images, sorted by reversal point.
savename = ['cp_aligned_images_newCsMinus_' expType];
ensureFigure(savename, 1);

params = struct();    
params.matpos = [0 0 1 1];    
params.figmargin = [0.15 0 0.1 0.2];
params.cellmargin = [0.025 0.025 0.025 0.025];    

hax = axesmatrix(3,2,1:6,params);

axes(hax(1)); hold on;
ylabel('Cue licks');
cData = good_licks;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsMinus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), markerSize, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

axes(hax(2)); hold on;
cData = aligned_licks(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', lineWidth);
set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

axes(hax(3)); hold on;
ylabel('ACh.', 'Color', mycolors('ChAT'));
cData = good_ch1;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsMinus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), markerSize, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

axes(hax(4)); hold on;
cData = aligned_ch1(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', lineWidth); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

axes(hax(5)); hold on;
ylabel('Dop.', 'Color', mycolors('DAT'));
cData = good_ch2;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsMinus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), markerSize, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', [0 10 20]);
xlabel('Trials from rev.');

axes(hax(6)); hold on;
cData = aligned_ch2(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', lineWidth); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', [-10 0 10]);
xlabel('Trials from c.p.');

formatFigurePublish('size', figSize, 'fontSize', 6);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
%     export_fig(fullfile(savepath, savename), '-eps');
end

%  averages
savename = ['cp_aligned_newCsMinus_avgs' '_' expType];
ensureFigure(savename, 1);
marker = '.';
yLim = [-0.8 0.8];

axes; hold on;


[hl, hp] = boundedline(xData, [nanmean(aligned_ch1)' nanmean(aligned_ch2)'], permute([nanSEM(aligned_ch1)' nanSEM(aligned_ch2)'], [1 3 2]),...
    'cmap', [mycolors('chat'); mycolors('dat')], 'nan', 'gap');
set(gca, 'YLim', yLim);
ylabel('Cue response (z-score)');
xlabel('Trials from lick c.p.');

formatFigurePublish('size', avgfigSize, 'fontSize', 6);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
%     export_fig(fullfile(savepath, savename), '-eps');
end



%% zscored averages
expType = 'all';
smoothWindow = 3;
compile_reversal_data;


trialWindow = [-30 30];
ylim = [-1 2];
figsize = [1.4 1.2];
% new cs plus
% common = sum(~isnan(newCsPlus.licks_cs)) > 3;
common = newCsPlus.firstRevTrial + trialWindow(1):newCsPlus.firstRevTrial + trialWindow(2) - 1;
savename = 'reversals_newCsPlus_zscored_all';
ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsPlus.trialNumber(common), nanmean(nanzscore(newCsPlus.licks_cs(goodReversals, common),0,2)), nanSEM(nanzscore(newCsPlus.licks_cs(goodReversals, common),0,2))',...
    'cmap', mycolors('licks'), 'nan', 'gap');
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus.trialNumber(common), nanmean(nanzscore(newCsPlus.phPeakMean_cs_ch2(goodReversals, common),0,2)), nanSEM(nanzscore(newCsPlus.phPeakMean_cs_ch2(goodReversals, common),0,2))',...
    'cmap', mycolors('dat'), 'nan', 'gap'); hold on
hla(2) = hl;
[hl, hp] = boundedline(newCsPlus.trialNumber(common), nanmean(nanzscore(newCsPlus.phPeakMean_cs_ch1(goodReversals, common),0,2)), nanSEM(nanzscore(newCsPlus.phPeakMean_cs_ch1(goodReversals, common),0,2))',...
    'cmap', mycolors('chat'), 'nan', 'gap');
hla(3) = hl;


set(hla, 'LineWidth', 1);
set(gca, 'XLim', trialWindow);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 1);
legend(hla, { '\bf\color[rgb]{0.5,0.5,0.5}Licks', '\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.'},...
    'Location', 'northwest', 'FontSize', 6, 'Interpreter', 'tex', 'Box', 'off');

title('New Cs+');
xlabel('Odor presentations from rev.');
ylabel('Cue response');

formatFigurePublish('size', figsize, 'fontSize', 6);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
end   


% new cs minus
% common = sum(~isnan(newCsMinus.licks_cs)) > 3;
common = newCsMinus.firstRevTrial + trialWindow(1):newCsMinus.firstRevTrial + trialWindow(2) - 1;

savename = 'reversals_newCsMinus_zscored_all';
ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsMinus.trialNumber(common), nanmean(nanzscore(newCsMinus.licks_cs(goodReversals, common),0,2)), nanSEM(nanzscore(newCsMinus.licks_cs(goodReversals, common),0,2))',...
    'cmap', mycolors('licks'));
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus.trialNumber(common), nanmean(nanzscore(newCsMinus.phPeakMean_cs_ch2(goodReversals, common),0,2)), nanSEM(nanzscore(newCsMinus.phPeakMean_cs_ch2(goodReversals, common),0,2))',...
    'cmap', mycolors('dat')); hold on
hla(2) = hl;
[hl, hp] = boundedline(newCsMinus.trialNumber(common), nanmean(nanzscore(newCsMinus.phPeakMean_cs_ch1(goodReversals, common),0,2)), nanSEM(nanzscore(newCsMinus.phPeakMean_cs_ch1(goodReversals, common),0,2))',...
    'cmap', mycolors('chat'));
hla(3) = hl;

set(hla, 'LineWidth', 1);
set(gca, 'XLim', trialWindow);%, 'YLim', ylim);    
h  = addOrginLines;
set(h, 'LineWidth', 1);
legend(hla, {'\bf\color[rgb]{0.5,0.5,0.5}Licks', '\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.'},...
            'Location', 'northeast', 'FontSize', 6, 'Interpreter', 'tex', 'Box', 'off');
title('New Cs-');
xlabel('Odor presentations from rev.');
ylabel('Cue response');    
formatFigurePublish('size', figsize, 'fontSize', 6);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));
end   

