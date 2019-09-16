%{
developed as subset to reversals_noPunish_poolReversals, focuses on trying
to leverage statistical power of paired recordings of ACh. and Dop. neurons

%}
DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, 'pooled', filesep);
smoothWindow = 1; % smoothing ruins cross correlation validity, also problematic for permutation test
saveOn = 1;

%%
exp.value = {'DC_44'  'DC_46'  'DC_47' 'DC_53'  'DC_54'  'DC_56'}; % exclude DC_51
exp.valence = {'DC_17'  'DC_20'  'DC_35'  'DC_36'  'DC_37'  'DC_40'};
exp.all = [exp.value exp.valence];
expType = 'all';


%%
compile_reversal_data;

%% changepoint detection for acquisition (acq) and extinction (ext)
% compFields = {'csPlus', 'csMinus'};
% fitFields = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'};

baselineTrials = 20;
% baselineTrials = 1;
% cp_licks = bpChangePoints([AR.csMinus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csPlus.licks_cs.after(:, 1:end)], 2, 1000);
% cp.csPlus.licks_cs = bpChangePoints(smoothdata(newCsPlus.licks_cs(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 'movmean', smoothWindow, 'omitnan'), 2, 1000, 'up');
cp.csPlus.licks_cs = bpChangePoints(newCsPlus.licks_cs(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
cp.csPlus.phPeakMean_cs_ch1 = bpChangePoints(newCsPlus.phPeakMean_cs_ch1(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
cp.csPlus.phPeakMean_cs_ch2 = bpChangePoints(newCsPlus.phPeakMean_cs_ch2(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
cp.csMinus.licks_cs = bpChangePoints(newCsMinus.licks_cs(:,newCsMinus.firstRevTrial - baselineTrials:end), 2, 1000, 'down');
cp.csMinus.phPeakMean_cs_ch1 = bpChangePoints(newCsMinus.phPeakMean_cs_ch1(:,newCsMinus.firstRevTrial - baselineTrials:end), 2, 1000, 'down');
cp.csMinus.phPeakMean_cs_ch2 = bpChangePoints(newCsMinus.phPeakMean_cs_ch2(:,newCsMinus.firstRevTrial - baselineTrials:end), 2, 1000, 'down');

% cp.csPlus.licks_cs = bpChangePoints([AR.csMinus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csPlus.licks_cs.after(:, 1:end)], 2, 1000, 'up');
% cp.csPlus.phPeakMean_cs_ch1 = bpChangePoints([AR.csMinus.phPeakMean_cs_ch1.before(:, end - baselineTrials + 1:end) AR.csPlus.phPeakMean_cs_ch1.after(:, 1:end)], 2, 1000, 'up');
% cp.csPlus.phPeakMean_cs_ch2 = bpChangePoints([AR.csMinus.phPeakMean_cs_ch2.before(:, end - baselineTrials + 1:end) AR.csPlus.phPeakMean_cs_ch2.after(:, 1:end)], 2, 1000, 'up');
% cp.csMinus.licks_cs = bpChangePoints([AR.csPlus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csMinus.licks_cs.after(:, 1:end)], 2, 1000, 'down');
% cp.csMinus.phPeakMean_cs_ch1 = bpChangePoints([AR.csPlus.phPeakMean_cs_ch1.before(:, end - baselineTrials + 1:end) AR.csMinus.phPeakMean_cs_ch1.after(:, 1:end)], 2, 1000, 'down');
% cp.csMinus.phPeakMean_cs_ch2 = bpChangePoints([AR.csPlus.phPeakMean_cs_ch2.before(:, end - baselineTrials + 1:end) AR.csMinus.phPeakMean_cs_ch2.after(:, 1:end)], 2, 1000, 'down');


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
%

goodPupil = auROC.pupil_csBaselined.before > 0.2;
goodWhisk = auROC.whisk_csBaselined.before > 0.2;


%% take advantage of paired recordings, subtract dopamine cue response from Ach. cue response, etc.
newCsPlus.phPeakMean_cs_AchMinusDop = newCsPlus.phPeakMean_cs_ch1 - newCsPlus.phPeakMean_cs_ch2;
newCsMinus.phPeakMean_cs_AchMinusDop = newCsMinus.phPeakMean_cs_ch1 - newCsMinus.phPeakMean_cs_ch2;
alwaysCsPlus.phPeakMean_cs_AchMinusDop = alwaysCsPlus.phPeakMean_cs_ch1 - alwaysCsPlus.phPeakMean_cs_ch2;
odor3.phPeakMean_cs_AchMinusDop = odor3.phPeakMean_cs_ch1 - odor3.phPeakMean_cs_ch2;


%% images
% orderings = {'newCsPlus', 'newCsMinus', 'alwaysCsPlus'};
fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'licks_cs', 'csLicksROC', 'pupil_cs', 'whisk_cs'};
titles = {'ACh', 'Dop.', 'Licks', 'Lick auROC', 'Pupil', 'Whisk'};
clim = [-5 5];
fh=[];

savename = ['newCsPlus_image_' expType];
    fh(end+1) = ensureFigure(savename, 1);
cLimFactor = 3;
% smoothWindow = 1;
xData = [min(newCsPlus.trialNumber), max(newCsPlus.trialNumber)];
xlim = [-30 70];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = newCsPlus.(sfield)(sortOrder, :);
%     cData = smoothdata(cData, 2, 'movmean', smoothWindow, 'omitnan');
    imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:length(sortOrder), [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 length(sortOrder)]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs+');
set(gcf, 'Position', [304   217   633   485]);

savename = ['newCsMinus_image_' expType];
fh(end+1) = ensureFigure(savename, 1);
cLimFactor = 3;
xData = [min(newCsMinus.trialNumber), max(newCsMinus.trialNumber)];
xlim = [-30 70];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = newCsMinus.(sfield)(sortOrder, :);
    cData = smoothdata(cData, 2, 'movmean', smoothWindow, 'omitnan');
    imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:length(sortOrder), [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 length(sortOrder)]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs-');
set(gcf, 'Position', [304   217   633   485]);



%% images, ch1, ch2, and licks ONLY, ordered by changepoints
% orderings = {'newCsPlus', 'newCsMinus', 'alwaysCsPlus'};
fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'licks_cs'};
titles = {'ACh', 'Dop.', 'Licks'};
clim = [-5 5];
fh=[];

savename = ['newCsPlus_image_ordered' '_' expType];
fh(end+1) = ensureFigure(savename, 1);
cLimFactor = 3;
xData = [min(newCsPlus.trialNumber), max(newCsPlus.trialNumber)];
xlim = [-30 70];
nFields = length(fieldsToShow);
for fcounter = 1:nFields
    sfield = fieldsToShow{fcounter};    

    for ocounter = 1:nFields
        ofield = fieldsToShow{ocounter};
        sortVariable = cp.csPlus.(ofield).index;
        sortVariable(~goodReversals) = NaN;
        [sorted, sortOrder] = sort(sortVariable);

        subplot(nFields,nFields,(ocounter - 1) * nFields + fcounter);
        cData = newCsPlus.(sfield)(sortOrder, :);
        cData = smoothdata(cData, 2, 'movmean', 1, 'omitnan');
        imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
        scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:length(sortOrder), [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
        set(gca, 'YLim', [1 length(sortOrder)]);
        set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
        t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
    end
end
subplot(nFields,nFields,nFields * nFields - floor(nFields/2)); xlabel('Odor presentations from reversal');
subplot(nFields,nFields,floor(nFields/2)); title('New Cs+');
set(gcf, 'Position', [304   217   633   485]);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

savename = ['newCsMinus_image_ordered' '_' expType];
fh(end+1) = ensureFigure(savename, 1);
cLimFactor = 3;
xData = [min(newCsMinus.trialNumber), max(newCsMinus.trialNumber)];
xlim = [-30 70];
nFields = length(fieldsToShow);
for fcounter = 1:nFields
    sfield = fieldsToShow{fcounter};    

    for ocounter = 1:nFields
        ofield = fieldsToShow{ocounter};
        sortVariable = cp.csMinus.(ofield).index;
        sortVariable(~goodReversals) = NaN;
        [sorted, sortOrder] = sort(sortVariable);

        subplot(nFields,nFields,(ocounter - 1) * nFields + fcounter);
        cData = newCsMinus.(sfield)(sortOrder, :);
        cData = smoothdata(cData, 2, 'movmean', 1, 'omitnan');
        imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
        scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:length(sortOrder), [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
        set(gca, 'YLim', [1 length(sortOrder)]);
        set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
        t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
    end
end
subplot(nFields,nFields,nFields * nFields - floor(nFields/2)); xlabel('Odor presentations from reversal');
subplot(nFields,nFields,floor(nFields/2)); title('New Cs-');
set(gcf, 'Position', [304   217   633   485]);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

%% averages
trialWindow = [-30 30];
ylim = [-1 2];
figsize = [1.7 1.5];
% new cs plus
% common = sum(~isnan(newCsPlus.licks_cs)) > 3;
common = newCsPlus.firstRevTrial + trialWindow(1):newCsPlus.firstRevTrial + trialWindow(2) - 1;
savename = ['reversals_newCsPlus' expType];
fh = ensureFigure(savename, 1);
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

formatFigurePublish('size', figsize);
if saveOn 
    export_fig(fullfile(savepath, savename), '-eps');
end   


% new cs minus
% common = sum(~isnan(newCsMinus.licks_cs)) > 3;
common = newCsMinus.firstRevTrial + trialWindow(1):newCsMinus.firstRevTrial + trialWindow(2) - 1;

savename = ['reversals_newCsMinus' expType];
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
formatFigurePublish('size', figsize);
if saveOn 
    export_fig(fullfile(savepath, savename), '-eps');
end   



%% test subtraction approach to take advantage of paired recordings

common = sum(~isnan(newCsPlus.licks_cs)) > 3;
common_odor3 = (-9 <= odor3_trialNumber) & (odor3_trialNumber <= 9);


savename = ['ACh_minus_Dop_avgs_' '_' expType];
ensureFigure(savename, 1);
trialRange = [-30 30];
subplot(2,2,1);
title('new Cs+');
[hl, hp] = boundedline(newCsPlus.trialNumber(common), nanmean(newCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common)), nanSEM(newCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common))',...
    'cmap', [0 1 0], 'nan', 'gap'); hold on
set(gca, 'XLim', trialRange);%, 'YLim', [-1 2]);
t = textBox('\color[rgb]{0,0.5,0}ACh. - Dop.'); 
set(t, 'Interpreter', 'tex', 'FontSize', 8, 'FontWeight', 'bold');
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');
addOrginLines;
% t = textBox('\color{green} ACh. minus Dop.'); set(t, 'Interpreter', 'tex', 'FontSize', 8);



common = sum(~isnan(newCsMinus.licks_cs)) > 3;
subplot(2,2,2);
title('new Cs-');
[hl, hp] = boundedline(newCsMinus.trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_AchMinusDop(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_AchMinusDop(goodReversals, common))',...
    'cmap', [0 1 0], 'nan', 'gap'); hold on
set(gca, 'XLim', trialRange);%, 'YLim', [-1 2]);
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');
addOrginLines;
% t = textBox('\color{green} ACh. minus Dop.'); set(t, 'Interpreter', 'tex', 'FontSize', 8);
subplot(2,2,3);
title('always Cs+');
[hl, hp] = boundedline(alwaysCsPlus.trialNumber(common), nanmean(alwaysCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common)), nanSEM(alwaysCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common))',...
    'cmap', [0 1 0], 'nan', 'gap'); hold on
set(gca, 'XLim', trialRange);%, 'YLim', [-1 2]);
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');
addOrginLines;
% t = textBox('\color{green} ACh. minus Dop.'); set(t, 'Interpreter', 'tex', 'FontSize', 8);

subplot(2,2,4);
title('odor3');
[hl, hp] = boundedline(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_AchMinusDop(goodReversals, common_odor3)), nanSEM(odor3.phPeakMean_cs_AchMinusDop(goodReversals, common_odor3))',...
    'cmap', [0 1 0], 'nan', 'gap'); hold on
% % set(gca, 'XLim', trialRange);%, 'YLim', [-1 2]);
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');
addOrginLines;

% t = textBox('\color[rgb]{0,0.5,0}ACh. - Dop.'); 
% set(t, 'Interpreter', 'tex', 'FontSize', 8, 'FontWeight', 'bold');
formatFigurePublish('size', [3 2]);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    


%% align reversals by lick changepoint
minLogit = 5;
trialWindow = [-20 20];
cpField = 'licks_cs';
cLimFactor = 4;
trialWindow_rev = [-10 30];
markerSize = 4;


% first new csPlus (acquisition)

% setup for aligned by changepoint
goodOnes = goodReversals & (cp.csPlus.(cpField).logit > minLogit) & (cp.csPlus.(cpField).index > (baselineTrials - 0));
zeroTrials = cp.csPlus.(cpField).index(goodOnes) - baselineTrials;
reversalPoints = 0 - zeroTrials;
reversalPoints = reversalPoints;
[sorted, ix] = sort(reversalPoints); % THEN sort them


% first select good subset

good_licks = newCsPlus.licks_cs(goodOnes, :);
% try normalizing licks
good_licks = good_licks ./ percentile(good_licks(:,newCsPlus.firstRevTrial:end), 0.9, 2);
good_ch1 = nanzscore(newCsPlus.phPeakMean_cs_ch1(goodOnes, :), 0, 2);
good_ch2 = nanzscore(newCsPlus.phPeakMean_cs_ch2(goodOnes, :), 0, 2);
good_subtract = newCsPlus.phPeakMean_cs_AchMinusDop(goodOnes, :);

reversals = true(sum(goodOnes), 1);
[aligned_subtract, xData] = alignedDataWindow(good_subtract, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1), sum(goodOnes), 1));
[aligned_licks, ~] = alignedDataWindow(good_licks, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch1, ~] = alignedDataWindow(good_ch1, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch2, ~] = alignedDataWindow(good_ch2, reversals, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));

% randomly permute across reversals (not time) to preserve average but
% scramble potential changepoints
maxTrials = length(newCsPlus.trialNumber);

% generate subscript indices for permutation
row_ix = repmat(1:maxTrials, sum(goodOnes), 1);
col_ix = zeros(sum(goodOnes), maxTrials);
for counter = 1:maxTrials
    col_ix(:,counter) = randperm(sum(goodOnes))';
end



lin_ix = sub2ind([sum(goodOnes) maxTrials], col_ix, row_ix);
perm_licks = good_licks(lin_ix);
perm_ch1 = good_ch1(lin_ix);
perm_ch2 = good_ch2(lin_ix);
perm_subtract = good_subtract(lin_ix);

% quick plot of regular and permuted lick_cs images, and corresponding
% averages
savename = 'permutation_sanity_check';
ensureFigure(savename, 1);
clim = [0 1]; xlim = [-20 40];
subplot(3,1,1);
imagesc(newCsPlus.trialNumber(:, newCsPlus.firstRevTrial - 21:newCsPlus.firstRevTrial + 39), 1:sum(goodOnes), good_licks(:, newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40), clim);
set(gca, 'XLim', xlim);
title('Licks');
subplot(3,1,2);
imagesc(newCsPlus.trialNumber(:, newCsPlus.firstRevTrial - 21:newCsPlus.firstRevTrial + 39), 1:sum(goodOnes), perm_licks(:, newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40), clim);
set(gca, 'XLim', xlim);
title('Licks permuted');
subplot(3,1,3); hold on; 
title('averages');
plot(newCsPlus.trialNumber(newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40), nanmean(good_licks(:, newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40)), '-k');
plot(newCsPlus.trialNumber(newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40), nanmean(perm_licks(:, newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40)), '--', 'Color', [0.8 0.8 0.8]);
legend({'intact', 'permuted'}, 'Location', 'best');
set(gca, 'XLim', xlim); xlabel('new Cs+ trials from rev.');

formatFigurePoster([3 5], '', 10);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
end    


% cp_perm.csPlus.licks_cs = bpChangePoints([AR.csMinus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csPlus.licks_cs.after(:, 1:end)], 2, 1000);
cp_perm.csPlus.licks_cs = bpChangePoints(perm_licks(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
zeroTrials_perm = cp_perm.csPlus.licks_cs.index - baselineTrials;
reversalPoints_perm = 0 - zeroTrials_perm;
reversalPoints_perm = reversalPoints_perm;
[sorted_perm, ix_perm] = sort(reversalPoints_perm); % THEN sort them

[aligned_subtract_perm, xData] = alignedDataWindow(perm_subtract, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1), sum(goodOnes), 1));
[aligned_licks_perm, ~] = alignedDataWindow(perm_licks, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch1_perm, ~] = alignedDataWindow(perm_ch1, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch2_perm, ~] = alignedDataWindow(perm_ch2, reversals, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));


% setup for alignment by reversal
% [sorted_rev, sortOrder_rev] = sort(cp.csPlus.licks_cs.index(goodOnes));
cp_rev = zeroTrials;
cp_rev_perm = zeroTrials_perm;

% first images, sorted by reversal point.
savename = ['cp_aligned_newCsPlus_images' '_' expType];
ensureFigure(savename, 1);

ckey = pink;
nColors = size(ckey, 1);

max_logit = max(cp.csPlus.licks_cs.logit(goodOnes & isfinite(cp.csPlus.licks_cs.logit)));
max_logit = max_logit * 1.2; % make Inf look brighter than others
color_ix = cp.csPlus.licks_cs.logit(goodOnes);
color_ix(~isfinite(color_ix)) = max_logit;
color_ix = ceil(color_ix ./ max_logit .* nColors);
cp_colors = ckey(color_ix, :);


color_ix_perm = cp_perm.csPlus.licks_cs.logit;
color_ix_perm(~isfinite(color_ix)) = max_logit;
color_ix_perm = ceil(color_ix ./ max_logit .* nColors);
cp_colors_perm = ckey(color_ix, :);



subplot(4,4,1); hold on;
title('Reversal Aligned');
ylabel('Cue licks');
cData = good_licks;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), 20, cp_colors, 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);

subplot(4,4,2); hold on;
title('Changepoint Aligned');
cData = aligned_licks(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4);
set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,3); hold on;
title('Rev. aligned, permuted');
cData = perm_licks;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev_perm, 1:sum(goodOnes), 20, cp_colors_perm, 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);


subplot(4,4,4); hold on;
title('CP aligned, permuted');
cData = aligned_licks_perm(ix_perm, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,5); hold on;
ylabel('ACh.', 'Color', mycolors('ChAT'));
cData = good_ch1;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), 20, cp_colors, 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

subplot(4,4,6); hold on;
cData = aligned_ch1(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,7); hold on;
cData = perm_ch1;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev_perm, 1:sum(goodOnes), 20, cp_colors_perm, 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);

subplot(4,4,8); hold on;
cData = aligned_ch1_perm(ix_perm, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,9); hold on;
ylabel('Dop.', 'Color', mycolors('DAT'));
cData = good_ch2;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), 20, cp_colors, 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

subplot(4,4,10); hold on;
cData = aligned_ch2(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,11); hold on;
cData = perm_ch2;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev_perm, 1:sum(goodOnes), 20, cp_colors_perm, 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);

subplot(4,4,12); hold on;
cData = aligned_ch2_perm(ix_perm, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,13); hold on;
ylabel('Subtract', 'Color', [0 1 0]);
cData = good_licks;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), 20, cp_colors, 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);
xlabel('Trials from reversal');

subplot(4,4,14); hold on;
cData = aligned_subtract(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', []);
xlabel('Trials from changepoint');

subplot(4,4,15); hold on;
cData = perm_subtract;
% use same clim for all images in each row
clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev_perm, 1:sum(goodOnes), 20, cp_colors_perm, 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);

subplot(4,4,16); hold on;
cData = aligned_subtract_perm(ix_perm, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', []);
xlabel('Trials from changepoint');



formatFigurePoster([8 10], '', 10);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

% second averages
savename = ['cp_aligned_newCsPlus_avgs' '_' expType];
ensureFigure(savename, 1);
marker = '.';
subplot(2,2,1); title('new Cs+ (acquisition)');
[hl, hp] = boundedline(xData, [nanmean(aligned_licks)' nanmean(aligned_licks_perm)'], permute([nanSEM(aligned_licks)' nanSEM(aligned_licks_perm)'], [1 3 2]),...
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

%
savename = 'cp_licks_cp_vs_permute_correlation';
ensureFigure(savename, 1);
title('Cs+ licking changepoints');
xdata = cp.csPlus.licks_cs.index(goodOnes) - baselineTrials;
ydata = cp_perm.csPlus.licks_cs.index - baselineTrials;
figure; scatter(xdata, ydata);
xlabel('change points'); ylabel('change points (permuted)');
textBox(sprintf('Rsq=%.3g', corr(xdata, ydata)));
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
end    

%% time to 50% or equivalent fraction

fract = 0.5;






%% weibull function
baselineTrials = 20;
fitField = 'phPeakMean_cs_ch1';


%% fit weibull function, CDF form with offset and scaling parameters

compFields = {'csPlus', 'csMinus';...  % first row is after reversal
              'csMinus', 'csPlus'...   % second row is before reversal
              };
fitFields = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'};
outputFields = {'object', 'gof', 'output', 'toFit', 'a', 'b', 'c', 'd'};
weibull = struct();

nReversals = size(AR.csPlus.globalTrialNumber.after, 1);

for counter = 1:length(compFields)
    for counter2 = 1:length(fitFields)
        for counter3 = 1:length(outputFields)
            if ~ismember(outputFields{counter3}, {'a', 'b', 'c', 'd'})
                weibull.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = cell(nReversals, 1);
            else
                weibull.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = NaN(nReversals, 1);
            end
        end
    end
end

weibullModel =  'a * (1 - exp(-1 * (x/b)^c)) + d'; % weibull function, CDF form


for compCounter = 1:size(compFields, 2)    
    for fieldCounter = 1:length(fitFields)
        weibullData = [AR.(compFields{2, compCounter}).(fitFields{fieldCounter}).before(:, end - baselineTrials + 1:end) AR.(compFields{1, compCounter}).(fitFields{fieldCounter}).after];
        for counter = 1:size(weibullData, 1)        
            toFit = weibullData(counter, ~isnan(weibullData(counter, :)));
            switch compFields{1, compCounter}
                case 'csPlus'
                    fo = fitoptions('Method', 'NonlinearLeastSquares',... 
                        'Upper', [Inf  Inf Inf Inf],...  % 20 (3rd upper)
                        'Lower', [0 0 0 -Inf],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
                        'StartPoint', [range(toFit) baselineTrials baselineTrials min(toFit)]...  % 'StartPoint', [range(toFit) baselineTrials baselineTrials min(toFit)]...
                        );
                case 'csMinus'
                    fo = fitoptions('Method', 'NonlinearLeastSquares',... 
                        'Upper', [0  Inf Inf Inf],...  % 20 (3rd upper)
                        'Lower', [-Inf 0 0 -Inf],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
                        'StartPoint', [-range(toFit) baselineTrials baselineTrials max(toFit)]... % 'StartPoint', [-range(toFit) baselineTrials baselineTrials max(toFit)]...
                        );            
            end            
            ft = fittype(weibullModel, 'options', fo);
            weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).toFit{counter} = toFit;
            try
                [fitobject, gof, output] = fit((0:length(toFit) - 1)', toFit', ft, fo);
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).object{counter} = fitobject;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).a(counter) = fitobject.a;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).b(counter) = fitobject.b;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).c(counter) = fitobject.c;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).d(counter) = fitobject.d;
                weibull.(compFields{1, compCounter}).(fitFields{fieldCounter}).output{counter} = output;
            catch
                continue
            end

        end
    end
end

%% fit exponentials to post-reversal data for both acquisition and extinction (no baseline trials for exponential....)
% ss = struct('object', zeros(sum(goodReversals), 1), 'gof', zeros(sum(goodReversals), 1), 'output', zeros(sum(goodReversals), 1), 'toFit', zeros(sum(goodReversals), 1));

compFields = {'csPlus', 'csMinus'};
fitFields = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'};
outputFields = {'object', 'gof', 'output', 'toFit', 'a', 'b', 'c'};
expFit = struct();
for counter = 1:length(compFields)
    for counter2 = 1:length(fitFields)
        for counter3 = 1:length(outputFields)
            if ~ismember(outputFields{counter3}, {'a', 'b', 'c'})
                expFit.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = cell(nReversals, 1);
            else
                expFit.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = NaN(nReversals, 1);
            end
        end
    end
end


% these options should work for newCsPlus and newCsMinus (up or down and
% for z scored and lick rate data (both should fall within interval of -100
% to 100
% fo = fitoptions('Method', 'NonlinearLeastSquares',...
%     'Upper', [100  100 1000],...
%     'Lower', [-100 -100 0],...    
%     'StartPoint', [0 1 10]...
%     );

fo = fitoptions('Method', 'NonlinearLeastSquares',...
    'Upper', [Inf  Inf Inf],...
    'Lower', [-Inf -Inf 0],...    
    'StartPoint', [0 1 10]...
    );

expModel = 'a + b * exp(-x/c)';
for compCounter = 1:length(compFields)
    for fieldCounter = 1:length(fitFields)
        expFitData = AR.(compFields{compCounter}).(fitFields{fieldCounter}).after;
        for counter = 1:size(expFitData, 1)        
            toFit = expFitData(counter, ~isnan(expFitData(counter, :)));
            
            switch compFields{compCounter}
                case 'csPlus'
                    fo = fitoptions('Method', 'NonlinearLeastSquares',...
                        'Upper', [Inf  0 Inf],...
                        'Lower', [-Inf -Inf 0],...    
                        'StartPoint', [0 -1 10]...
                        );
                case 'csMinus'
                    fo = fitoptions('Method', 'NonlinearLeastSquares',...
                        'Upper', [Inf  Inf Inf],...
                        'Lower', [-Inf 0 0],...    
                        'StartPoint', [0 1 10]...
                        );
            end                        
            ft = fittype(expModel, 'options', fo);
            expFit.(compFields{compCounter}).(fitFields{fieldCounter}).toFit{counter} = toFit;
            try
                [fitobject, gof, output] = fit((0:length(toFit) - 1)', toFit', ft, fo);
                expFit.(compFields{compCounter}).(fitFields{fieldCounter}).object{counter} = fitobject;
                expFit.(compFields{compCounter}).(fitFields{fieldCounter}).a(counter) = fitobject.a;
                expFit.(compFields{compCounter}).(fitFields{fieldCounter}).b(counter) = fitobject.b;
                expFit.(compFields{compCounter}).(fitFields{fieldCounter}).c(counter) = fitobject.c;
            catch
                continue
            end

        end
    end
end





%% make a bar graph or box plot of changepoints
    
savename = ['changepoints_all' '_' expType];
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


all_cps = [];

for counter = 1:size(fields, 1)
    ydata = cp.(fields{counter, 1}).(fields{counter, 2}).index - baselineTrials;       
    all_cps = [all_cps ydata];
%     scatter(repmat(counter, numel(ydata), 1) + (rand(numel(ydata), 1) - 0.5)/3, ydata, markerSize, fields{counter, 3}, '.');    
%     errorbar(counter, mean(ydata), std(ydata)/sqrt(numel(ydata)), 'Color', fields{counter, 3}, 'LineWidth', 2)
end
bar(1:6, median(all_cps), 'CData', vertcat(fields{:,3}), 'EdgeColor', 'flat', 'FaceColor', 'none'); 
violins = violinplot(all_cps, fields(:, 4), 'ShowMean', false, 'ShowNotches', false, 'EdgeColor', [1 1 1]);
for counter = 1:length(violins)
    violins(counter).ViolinColor = fields{counter,3};
    xData = violins(counter).ScatterPlot.XData;
    yData = violins(counter).ScatterPlot.YData;
    violins(counter).ScatterPlot.delete;
    violins(counter).ScatterPlot = scatter(xData, yData, 12, fields{counter, 3}, 'filled');    
    violins(counter).MedianPlot.delete;
end
ylabel('changepoint trials from rev.');
% set(gca, 'XLim', [0.5 6.5], 'XTick', 1:6, 'XTickLabel', fields(:,4)', 'FontSize', 12); ylabel('trials from rev.', 'FontSize', 12);
% set(gca, 'YLim', [0 20]);
% 
% for counter = 1:size(fields, 1)
%     ydata = cp.(fields{counter, 1}).(fields{counter, 2}).index(goodReversals) - baselineTrials;   
%     kstest(ydata)
%     all_cps = [all_cps ydata];
%     scatter(repmat(counter, numel(ydata), 1) + (rand(numel(ydata), 1) - 0.5)/3, ydata, markerSize, fields{counter, 3}, '.');
%     errorbar(counter, mean(ydata), std(ydata)/sqrt(numel(ydata)), 'Color', fields{counter, 3}, 'LineWidth', 2)
% end
% 
% set(gca, 'XLim', [0.5 6.5], 'Ylim', [-20 40], 'XTick', 1:6, 'XTickLabel', fields(:,4)', 'FontSize', 12); ylabel('trials from rev.', 'FontSize', 12);
plot(get(gca, 'XLim'), [0 0], '--', 'Color', [0.7 0.7 0.7]);
set(gca, 'YLim', [-baselineTrials 40]);
formatFigurePublish('size', [2.6 1.5], 'fontSize', 8);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    export_fig(fullfile(savepath, savename), '-eps');
end
%% stats on changepoints


comp = {'Lick_vs_ACh'; 'Lick_vs_Dop'; 'ACh_vs_Dop'};
acq_p = zeros(3,1);
ext_p = zeros(3,1);
cp_stats = table(comp, acq_p, ext_p);
cp_stats.acq_p(1) = signrank(all_cps(goodReversals,1), all_cps(goodReversals,2));
cp_stats.acq_p(2) = signrank(all_cps(goodReversals,1), all_cps(goodReversals,3));
cp_stats.acq_p(3) = signrank(all_cps(goodReversals,2), all_cps(goodReversals,3));
cp_stats.ext_p(1) = signrank(all_cps(goodReversals,4), all_cps(goodReversals,5));
cp_stats.ext_p(2) = signrank(all_cps(goodReversals,4), all_cps(goodReversals,6));
cp_stats.ext_p(3) = signrank(all_cps(goodReversals,5), all_cps(goodReversals,6));

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



%% plot example reversals with weibull fits, changepoints
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
    
    % weibull
    subplot(nShow, 3, counter*3 - 2); hold on;
    if counter == 1
        title('Weibull');
    end
    plot(weibull.(condition).(dataField).toFit{thisRev}, '.', 'Color', [0 0.75 0]); 
    plot(weibull.(condition).(dataField).object{thisRev}); legend off;
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
    scatter(1:length(cp.(condition).(dataField).cumsum{thisRev}), cp.(condition).(dataField).cumsum{thisRev}, 10, cp.(condition).(dataField).logitAll{thisRev}); colormap jet;
    line(repmat(cp.(condition).(dataField).index(thisRev), 1, 2), get(gca, 'YLim')); 
    plot([1 1] + baselineTrials, get(gca, 'YLim'), '--k');
    set(gca, 'XLim', [0 120]);
    textBox(sprintf('Logit=%.2f', cp.(condition).(dataField).logit(thisRev)));
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
% sortVariable = (cp_licks.index - baselineTrials); % e.g. 20 baseline trials
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



%% cross correlations, all mouse/odor/reversal pairs

savename = 'xcorr_all_corrected';
ensureFigure(savename, 1); 
% [R, lags] = avgXCorr(x, y, maxlag)
maxlag = 10;
trialRange = [-30 30];
R = deal(zeros(maxlag*2+1, 6, 3)); % corrected, shift predictor, raw occupy third dimension
lags = -maxlag:maxlag;

%     chat = nanzscore(newCsPlus.phPeakMean_cs_ch1, 0, 2);
%     dat = nanzscore(newCsPlus.phPeakMean_cs_ch2, 0, 2);
%     licks = nanzscore(newCsPlus.licks_cs, 0, 2);

chat = [newCsPlus.phPeakMean_cs_ch1(goodReversals, newCsPlus.firstRevTrial + trialRange(1):newCsPlus.firstRevTrial + trialRange(2) - 1);...
    newCsMinus.phPeakMean_cs_ch1(goodReversals, newCsMinus.firstRevTrial + trialRange(1):newCsMinus.firstRevTrial + trialRange(2) - 1)];
dat = [newCsPlus.phPeakMean_cs_ch2(goodReversals, newCsPlus.firstRevTrial + trialRange(1):newCsPlus.firstRevTrial + trialRange(2) - 1);...
    newCsMinus.phPeakMean_cs_ch2(goodReversals, newCsMinus.firstRevTrial + trialRange(1):newCsMinus.firstRevTrial + trialRange(2) - 1)];
licks = [newCsPlus.licks_cs(goodReversals, newCsPlus.firstRevTrial + trialRange(1):newCsPlus.firstRevTrial + trialRange(2) - 1);...
    newCsMinus.licks_cs(goodReversals, newCsMinus.firstRevTrial + trialRange(1):newCsMinus.firstRevTrial + trialRange(2) - 1)];        

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


% if saveOn
%     saveas(gcf, fullfile(savepath, [savename '.fig']));
%     saveas(gcf, fullfile(savepath, [savename '.jpg']));   
%     saveas(gcf, fullfile(savepath, [savename '.epsc']));   
% end

%% cross correlations, always cs+

savename = 'xcorr_alwaysCsPlus_corrected';
ensureFigure(savename, 1); 
% [R, lags] = avgXCorr(x, y, maxlag)
maxlag = 10;
trialRange = [30 30];
R = deal(zeros(maxlag*2+1, 6, 3)); % corrected, shift predictor, raw occupy third dimension
lags = -maxlag:maxlag;

%     chat = nanzscore(newCsPlus.phPeakMean_cs_ch1, 0, 2);
%     dat = nanzscore(newCsPlus.phPeakMean_cs_ch2, 0, 2);
%     licks = nanzscore(newCsPlus.licks_cs, 0, 2);

chat = [alwaysCsPlus.phPeakMean_cs_ch1(goodReversals, alwaysCsPlus.firstRevTrial - trialRange:alwaysCsPlus.firstRevTrial + trialRange - 1);...
    newCsMinus.phPeakMean_cs_ch1(goodReversals, newCsMinus.firstRevTrial - trialRange:newCsMinus.firstRevTrial + trialRange - 1)];
dat = [alwaysCsPlus.phPeakMean_cs_ch2(goodReversals, alwaysCsPlus.firstRevTrial - trialRange:alwaysCsPlus.firstRevTrial + trialRange - 1);...
    newCsMinus.phPeakMean_cs_ch2(goodReversals, newCsMinus.firstRevTrial - trialRange:newCsMinus.firstRevTrial + trialRange - 1)];
licks = [alwaysCsPlus.licks_cs(goodReversals, alwaysCsPlus.firstRevTrial - trialRange:alwaysCsPlus.firstRevTrial + trialRange - 1);...
    newCsMinus.licks_cs(goodReversals, newCsMinus.firstRevTrial - trialRange:newCsMinus.firstRevTrial + trialRange - 1)];       

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
    
   

subplot(3,2,1); plot(lags, squeeze(R(:,1,:))); title('chat vs dat'); %legend('corrected', 'shift predictor', 'raw'); legend boxoff;
subplot(3,2,2); plot(lags, squeeze(R(:,2,:))); title('chat');  %legend('corrected', 'shift predictor', 'raw'); legend boxoff;
subplot(3,2,3); plot(lags, squeeze(R(:,3,:))); title('dat'); 
subplot(3,2,4); plot(lags, squeeze(R(:,4,:))); title('chat vs licks');
subplot(3,2,5); plot(lags, squeeze(R(:,5,:))); title('dat vs licks'); xlabel('always cs+ trials');
subplot(3,2,6); plot(lags, squeeze(R(:,6,:))); title('licks'); xlabel('always cs+ trials');


% if saveOn
%     saveas(gcf, fullfile(savepath, [savename '.fig']));
%     saveas(gcf, fullfile(savepath, [savename '.jpg']));   
%     saveas(gcf, fullfile(savepath, [savename '.epsc']));   
% end


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
linecolors = repmat(linecolors, 10,1);
keepers = [1:6];

nMice = max(mouseNumber);
savename = 'reversal_trial_correlations';
ensureFigure(savename, 1);
% allTrials = true; 
ll = [];
subplot(1,2,1); hold on; subplot(1,2,2); hold on;
R_cs = zeros(nMice, 1);
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
    R_cs(counter) = corr(xData, yData);
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

% subplot(1,2,1); legend(ll, DB.animals(keepers), 'Interpreter', 'none');
subplot(1,2,1); legend(ll, DB.animals, 'Interpreter', 'none');
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
scatter(1:length(cp.csPlus.licks_cs.cumsum{thisRev}), cp.csPlus.licks_cs.cumsum{thisRev}, 10, cp.csPlus.licks_cs.logitAll{thisRev}); colormap jet;

line(repmat(cp.csPlus.licks_cs.index(thisRev), 1, 2), get(gca, 'YLim'), 'Color', [0 1 0]); 
plot([1 1] + baselineTrials, get(gca, 'YLim'), '--k');
set(gca, 'XLim', [0 120]);
t = textBox(sprintf('Logit=%.2f', cp.csPlus.licks_cs.logit(thisRev)));
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
% zeroTrials = cp.csPlus.licks_cs.index;
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

