DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);
saveOn = 1;
smoothWindow = 1;

%% create structure containing all reversals-  SKIPS DC_51
AR = struct('csPlus', [], 'csMinus', [], 'csPlusReward', [], 'thirdOdor', []); % AR = all reversals
for si = 1:length(DB.animals)    
    animal = DB.animals{si};
    if strcmp(animal, 'DC_51')
        continue;
    end
    load(fullfile(DB.path, 'pooled', ['RE_' animal '.mat']));
    for group = fieldnames(AR)'
        sgroup = group{:};
        for field = fieldnames(RE.(sgroup))'
            sfield = field{:};
            if any(strcmp(sfield, {'trialsBefore', 'trialsAfter'})) % these are special fields that don't contain before and after data
                continue
            end
            if si == 1
                AR.(sgroup).(sfield).before = RE.(sgroup).(sfield).before;
                AR.(sgroup).(sfield).after = RE.(sgroup).(sfield).after;
            else
                AR.(sgroup).(sfield).before = expandVertCat(AR.(sgroup).(sfield).before, RE.(sgroup).(sfield).before, 'right');
                AR.(sgroup).(sfield).after = expandVertCat(AR.(sgroup).(sfield).after, RE.(sgroup).(sfield).after, 'left');
            end
        end
    end
end

firstReversals = cellfun(@(x,y) ~strcmp(x(1:5), y(1:5)), AR.csPlus.filename.before(1:end-1,end), AR.csPlus.filename.before(2:end,end));
firstReversals = [false; firstReversals];
nReversals = size(AR.csPlus.globalTrialNumber.after, 1);
% reversal #s
revNumber = ones(nReversals,1);
mouseNumber = ones(nReversals,1);

thisRev = 1;
thisMouse = 1;
for counter = 2:nReversals
    if strcmp(AR.csPlus.filename.before{counter - 1,end}(1:5), AR.csPlus.filename.before{counter,end}(1:5))
        thisRev = thisRev + 1;
    else
        thisRev = 1;
        thisMouse = thisMouse + 1;
    end
    revNumber(counter) = thisRev;
    mouseNumber(counter) = thisMouse;
end
sortVariable = revNumber;
% sortVariable = auROC.phPeakMean_cs_ch1.after;
[~, sortOrder] = sort(sortVariable);





% quality control- calculate auROC and dPrime for relevent comparisons
trialWindow = [-20 60];%?

% initialize
comp = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'pupil_csBaselined', 'whisk_csBaselined'}; % comparisons
all_ways = struct(...
    'before', zeros(nReversals, 1),... % compare cs- and cs+ before reversal
    'after', zeros(nReversals, 1),... % compare cs- and cs+ after reversal
    'acq', zeros(nReversals, 1),... % acquisition: compare cs- before and cs+ after reversal
    'ext', zeros(nReversals, 1)...  % extinction: compare cs+ before and cs- after reversal
    );
% clear dPrime auROC
for field = comp
    dPrime.(field{:}) = all_ways;
    auROC.(field{:}) = all_ways;
end

% calculate metrics
for field = comp
    for rev = 1:nReversals
        % before
        dPrime.(field{:}).before(rev) = dPrime_SNR(AR.csPlus.(field{:}).before(rev,trialWindow(1) + end + 1:end), AR.csMinus.(field{:}).before(rev,trialWindow(1) + end + 1:end));
        auROC.(field{:}).before(rev) = rocarea(stripNaNs(AR.csPlus.(field{:}).before(rev,trialWindow(1) + end + 1:end)), stripNaNs(AR.csMinus.(field{:}).before(rev,trialWindow(1) + end + 1:end)), 'scale');
        % after
        dPrime.(field{:}).after(rev) = dPrime_SNR(AR.csPlus.(field{:}).after(rev,1:trialWindow(2)), AR.csMinus.(field{:}).after(rev,1:trialWindow(2)));
        auROC.(field{:}).after(rev) = rocarea(stripNaNs(AR.csPlus.(field{:}).after(rev,1:trialWindow(2))), stripNaNs(AR.csMinus.(field{:}).after(rev,1:trialWindow(2))), 'scale');        
        % acquisition
        dPrime.(field{:}).acq(rev) = dPrime_SNR(AR.csPlus.(field{:}).after(rev,1:trialWindow(2)), AR.csMinus.(field{:}).before(rev,trialWindow(1) + end + 1:end));
        auROC.(field{:}).acq(rev) = rocarea(stripNaNs(AR.csPlus.(field{:}).after(rev,1:trialWindow(2))), stripNaNs(AR.csMinus.(field{:}).before(rev,trialWindow(1) + end + 1:end)), 'scale');                
        % extinction
        dPrime.(field{:}).ext(rev) = dPrime_SNR(AR.csPlus.(field{:}).before(rev,trialWindow(1) + end + 1:end), AR.csMinus.(field{:}).after(rev,1:trialWindow(2)));
        auROC.(field{:}).ext(rev) = rocarea(stripNaNs(AR.csPlus.(field{:}).before(rev,trialWindow(1) + end + 1:end)), stripNaNs(AR.csMinus.(field{:}).after(rev,1:trialWindow(2))), 'scale');
    end    
end

% quality control % 2
% detect when auROC values first exceed threshold
rocThresh = 0.5;
trialsToCriterion = NaN(nReversals, 1);
% looping is just easier
for counter = 1:nReversals    
    thisRev = AR.csPlus.csLicksROC.after(counter, :);
    thisRev = thisRev > rocThresh;
    nt = find(thisRev, 1);
    if ~isempty(nt)
        trialsToCriterion(counter) = nt;
    end
end



%% filter reversals according to quality
goodReversals = ...
    ~isnan(trialsToCriterion) &...
    auROC.phPeakMean_cs_ch1.acq > 0 &...
    auROC.phPeakMean_cs_ch2.acq > 0;
sortVariable = trialsToCriterion;
sortVariable(~goodReversals) = NaN;



[sorted, sortOrder] = sort(sortVariable);
% sortOrder = sortOrder(~isnan(sorted));
%



goodPupil = auROC.pupil_csBaselined.before > 0.2;
goodWhisk = auROC.whisk_csBaselined.before > 0.2;
%% compile data
fieldsToCompile = {...
    'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'phPeakPercentile_cs_ch1', 'phPeakPercentile_cs_ch2', 'csLicksROC', 'licks_cs', 'pupil_cs', 'pupil_csBaselined' 'whisk_cs', 'wheel_baseline',...
    'phPeakMean_us_ch1_deconv', 'phPeakMean_us_ch2_deconv', 'phPeakPercentile_us_ch1_deconv', 'phPeakPercentile_us_ch2_deconv'};

newCsPlus = struct();
newCsMinus = struct();
alwaysCsPlus = struct();
alwaysCsPlusReward = struct();
odor3 = struct();

for counter = 1:length(fieldsToCompile)
    field = fieldsToCompile{counter};
    newCsPlus.(field) = smoothdata([AR.csMinus.(field).before AR.csPlus.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    newCsMinus.(field) = smoothdata([AR.csPlus.(field).before AR.csMinus.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    alwaysCsPlus.(field) = smoothdata([AR.csPlus.(field).before AR.csPlus.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    alwaysCsPlusReward.(field) = smoothdata([AR.csPlusReward.(field).before AR.csPlusReward.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    odor3.(field) = smoothdata([AR.thirdOdor.(field).before AR.thirdOdor.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
end

newCsPlus_trialNumber = (1:size(newCsPlus.licks_cs, 2)) - size(AR.csMinus.licks_cs.before, 2);
newCsPlus_firstRevTrial = size(AR.csMinus.licks_cs.before, 2) + 1;
newCsMinus_trialNumber = (1:size(newCsMinus.licks_cs, 2)) - size(AR.csPlus.licks_cs.before, 2);
newCsMinus_firstRevTrial = size(AR.csPlus.licks_cs.before, 2) + 1;
alwaysCsPlus_trialNumber = (1:size(alwaysCsPlus.licks_cs, 2)) - size(AR.csPlus.licks_cs.before, 2);
alwaysCsPlus_firstRevTrial = size(AR.csPlus.licks_cs.before, 2) + 1;
alwaysCsPlusReward_trialNumber = (1:size(alwaysCsPlusReward.licks_cs, 2)) - size(AR.csPlusReward.licks_cs.before, 2);
alwaysCsPlusReward_firstRevTrial = size(AR.csPlusReward.licks_cs.before, 2) + 1;
odor3_trialNumber = (1:size(odor3.licks_cs, 2)) - size(AR.thirdOdor.licks_cs.before, 2);

oldCsPlus_trialNumber = -(size(AR.csPlus.licks_cs.before, 2) - 1) : 0;
oldCsMinus_trialNumber = -(size(AR.csMinus.licks_cs.before, 2) - 1) : 0;

%% images
% orderings = {'newCsPlus', 'newCsMinus', 'alwaysCsPlus'};
fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'licks_cs', 'csLicksROC', 'pupil_cs', 'whisk_cs'};
titles = {'ACh', 'Dop.', 'Licks', 'Lick auROC', 'Pupil', 'Whisk'};
clim = [-5 5];
fh=[];
imageArraySize = [4 3];
saveName = 'newCsPlus_image';
h = ensureFigure(saveName, 1);
cLimFactor = 3;
xData = [min(newCsPlus_trialNumber), max(newCsPlus_trialNumber)];
xlim = [-30 70];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = newCsPlus.(sfield)(sortOrder, :);
    cData = smoothdata(cData, 2, 'movmean', 5, 'omitnan');
    imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:length(sortOrder), [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 length(sortOrder)]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 6, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs+');

formatFigurePublish('size', imageArraySize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

saveName = 'newCsMinus_image';
h = ensureFigure(saveName, 1);
cLimFactor = 3;
xData = [min(newCsMinus_trialNumber), max(newCsMinus_trialNumber)];
xlim = [-30 70];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = newCsMinus.(sfield)(sortOrder, :);
    cData = smoothdata(cData, 2, 'movmean', 5, 'omitnan');
    imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:length(sortOrder), [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 length(sortOrder)]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 6, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs-');

formatFigurePublish('size', imageArraySize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% to detect latency to learning, fit weibull function, also changepoint detection

% weibull doesn't work great due to overfitting and noisy data
baselineTrials = 20;

fitField = 'licks_cs';
model = 'a * (1 - exp(-1 * (x/b)^c)) + d'; % weibull function, CDF form
fitData = [AR.csMinus.(fitField).before(:, end - baselineTrials + 1:end) AR.csPlus.(fitField).after(:, :)];% - nanmean(AR.csMinus.licks_cs.before(goodReversals, end - baselineTrials + 1:end), 2); % zero/baseline data at start

weibull = struct('object', [], 'gof', [], 'output', [], 'toFit', [], 'kappa', []);
weibull = repmat(weibull, size(fitData, 1), 1);
for counter = 1:size(fitData, 1)
    toFit = fitData(counter, ~isnan(fitData(counter, :)));
    fo = fitoptions('Method', 'NonlinearLeastSquares',...
        'Upper', [Inf  Inf 20 Inf],...
        'Lower', [0 0 0 -Inf],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
        'StartPoint', [mean(toFit) baselineTrials baselineTrials min(toFit)]...
        );
    ft = fittype(model, 'options', fo);
%     xData = (0:length(toFit) - 1)';
    [fitobject, gof, output] = fit((0:length(toFit) - 1)', toFit', ft, fo);
    weibull(counter).object = fitobject;
    weibull(counter).gof = gof;
    weibull(counter).output = output;
    weibull(counter).toFit = toFit;
    coeffs = coeffvalues(fitobject);
    weibull(counter).kappa = coeffs(2);
%     weibull(counter).xData = xData - baselineTrials;
end

% changepoint detection
baselineTrials = 20;
cp_licks = bpChangePoints([AR.csMinus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csPlus.licks_cs.after(:, 1:end)], 2, 1000);

%% plot example reversals with weibull fits, changepoints
nShow = 4;
% good revs include [77 3 48 74];,   3 and 74 best examples
% mediocre,  84 and 91 and 42 are mediocre ones
% 4 is good but gradual
toShow = find(goodReversals);
toShow = toShow(randperm(sum(goodReversals), nShow));

%%
% nShow = 4;  toShow = [74 4 91 42];
saveName = 'weibull_vs_changepoint_dev';
ensureFigure(saveName, 1);

for counter = 1:nShow   
    thisRev = toShow(counter);
    % weibull
    subplot(nShow, 2, counter*2 - 1);

    plot(weibull(thisRev).toFit, 'g.'); hold on;
    plot(weibull(thisRev).object, 'predfunc'); legend off;
    set(gca, 'XLim', [0 120], 'XTick', [20 70], 'XTickLabel', {'0', '50'}); ylabel('antic. licks (1/s)'); xlabel('Cs+ trials from rev.');
    if counter == 1
        title('Weibull fit to Cs+ licking');
    end        
    % changepoint
    subplot(nShow, 2, counter * 2); hold on;
    if counter == 1
        title('Changepoint, logit by permutation');
    end    
    scatter(1:length(cp_licks.cumsum{thisRev}), cp_licks.cumsum{thisRev}, 10, cp_licks.logitAll{thisRev}); colormap jet;
    line(repmat(cp_licks.index(thisRev), 1, 2), get(gca, 'YLim')); 
    textBox(sprintf('Logit=%.2f', cp_licks.logit(thisRev)), [], [.7 0.4], 6);
    set(gca, 'XLim', [0 120], 'XTick', [20 70], 'XTickLabel', {'0', '50'}); ylabel('cum. lick rate'); xlabel('Cs+ trials from rev.');
end

formatFigurePublish('size', [4 4]);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end


%%
saveName = 'weibull_vs_changepoint_scatter';
ensureFigure(saveName, 1); colormap('winter');
line([0 20 20], [20 20 0], 'Color', 'r'); hold on;
scatter([weibull(:).kappa], cp_licks.index(:), 10, repmat([0.7 0.7 0.7], size(fitData, 1), 1) .* repmat(~goodReversals, 1, 3));
xlabel('Kappa (weibull)'); ylabel('change point');
set(gca, 'YLim', [0 70], 'XLim', [0 70]);
addUnityLine;
legend('reversal trial');
formatFigurePublish('size', [2 2]);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
    export_fig(fullfile(savepath, saveName), '-jpg');
end

%% images ordered by lick changepoints

fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'licks_cs', 'csLicksROC', 'pupil_cs', 'whisk_cs'};
titles = {'ACh', 'Dop.', 'Licks', 'Lick auROC', 'Pupil', 'Whisk'};
clim = [-5 5];
fh=[];

sortVariable = (cp_licks.index - baselineTrials); % e.g. 20 baseline trials
sortVariable(~goodReversals) = NaN;
sortVariable(cp_licks.logit <= 3) = NaN;
[sorted, sortOrder] = sort(sortVariable);


saveName = 'newCsPlus_image_changepoints';
ensureFigure(saveName, 1); colormap parula;
cLimFactor = 3;
xData = [min(newCsPlus_trialNumber), max(newCsPlus_trialNumber)];
xlim = [-30 70];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = newCsPlus.(sfield)(sortOrder, :);
    cData = smoothdata(cData, 2, 'movmean', 5, 'omitnan');
    imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    line(sorted, 1:length(sortOrder), 'Color', 'w', 'LineWidth', 2); 
    set(gca, 'YLim', [1 find(isnan(sorted), 1) - 1]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 6, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs+, ordered by changepoint (logit > 3)');
% imageArraySize
formatFigurePublish('size', imageArraySize);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

%%
% averages
common_odor3 = (-9 <= odor3_trialNumber) & (odor3_trialNumber <= 9);
xlim = [-30 30];
ylim = [-1 2];
figsize = [1.7 1.5];
% new cs plus
common = sum(~isnan(newCsPlus.licks_cs)) > 3;
savename = 'reversals_newCsPlus';
fh(end + 1) = ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus.licks_cs(goodReversals, common)), nanSEM(newCsPlus.licks_cs(goodReversals, common))',...
    'cmap', mycolors('licks'), 'nan', 'gap');
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus.phPeakMean_cs_ch2(goodReversals, common)), nanSEM(newCsPlus.phPeakMean_cs_ch2(goodReversals, common))',...
    'cmap', mycolors('dat'), 'nan', 'gap'); hold on
hla(2) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus.phPeakMean_cs_ch1(goodReversals, common)), nanSEM(newCsPlus.phPeakMean_cs_ch1(goodReversals, common))',...
    'cmap', mycolors('chat'), 'nan', 'gap');
hla(3) = hl;


set(hla, 'LineWidth', 1);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
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
common = sum(~isnan(newCsMinus.licks_cs)) > 3;

savename = 'reversals_newCsMinus';
fh(end + 1) = ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.licks_cs(goodReversals, common)), nanSEM(newCsMinus.licks_cs(goodReversals, common))',...
    'cmap', mycolors('licks'));
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_ch2(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_ch2(goodReversals, common))',...
    'cmap', mycolors('dat')); hold on
hla(2) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_ch1(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_ch1(goodReversals, common))',...
    'cmap', mycolors('chat'));
hla(3) = hl;

set(hla, 'LineWidth', 1);
set(gca, 'XLim', xlim);%, 'YLim', ylim);    
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







