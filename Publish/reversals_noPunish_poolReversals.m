
DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, 'pooled', filesep);
smoothWindow = 1;
saveOn = 1;
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
    'phPeakMean_us_ch1_deconv', 'phPeakMean_us_ch2_deconv', 'phPeakPercentile_us_ch1_deconv', 'phPeakPercentile_us_ch2_deconv', 'phBaseline_ch1', 'phBaseline_ch2'};

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

%% plot quality control metrics
colorVar = revNumber;

% auROC vs dPrime
for field = comp
    field = field{:};
    savename = ['auROC_vs_dPrime_Scatter_' field];
    ensureFigure(savename, 1); colormap jet
    subplot(2,2,1); scatter(auROC.(field).before, dPrime.(field).before, [], colorVar); ylabel('dPrime'); title('before'); addOrginLines;
    subplot(2,2,2); scatter(auROC.(field).after, dPrime.(field).after, [], colorVar); title('after'); addOrginLines;
    subplot(2,2,3); scatter(auROC.(field).acq, dPrime.(field).acq, [], colorVar); ylabel('dPrime'); xlabel('auROC'); title('acquisition'); addOrginLines;
    subplot(2,2,4); scatter(auROC.(field).ext, dPrime.(field).ext, [], colorVar); xlabel('auROC'); title('extinction'); addOrginLines;
end

% revNumber vs auROC
for field = comp
    field = field{:};
    savename = ['revNumber_vs_auROC_Scatter_' field];
    ensureFigure(savename, 1); colormap jet
    subplot(2,2,1); scatter(revNumber, auROC.(field).before, [], colorVar); ylabel('auROC'); title('before'); addOrginLines;
    subplot(2,2,2); scatter(revNumber, auROC.(field).after, [], colorVar); title('after'); addOrginLines;
    subplot(2,2,3); scatter(revNumber, auROC.(field).acq, [], colorVar); ylabel('auROC'); xlabel('rev #'); title('acquisition'); addOrginLines;
    subplot(2,2,4); scatter(revNumber, auROC.(field).ext, [], colorVar); xlabel('rev #'); title('extinction'); addOrginLines;
end

% mouseNumber vs auROC
for field = comp
    field = field{:};
    savename = ['mouseNumber_vs_auROC_Scatter_' field];
    ensureFigure(savename, 1); colormap jet
    subplot(2,2,1); scatter(mouseNumber, auROC.(field).before, [], colorVar); ylabel('auROC'); title('before'); addOrginLines;
    subplot(2,2,2); scatter(mouseNumber, auROC.(field).after, [], colorVar); title('after'); addOrginLines;
    subplot(2,2,3); scatter(mouseNumber, auROC.(field).acq, [], colorVar); ylabel('auROC'); xlabel('mouse #'); title('acquisition'); addOrginLines;
    subplot(2,2,4); scatter(mouseNumber, auROC.(field).ext, [], colorVar); xlabel('mouse #'); title('extinction'); addOrginLines;
end
%% images
% orderings = {'newCsPlus', 'newCsMinus', 'alwaysCsPlus'};
fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'licks_cs', 'csLicksROC', 'pupil_cs', 'whisk_cs'};
titles = {'ACh', 'Dop.', 'Licks', 'Lick auROC', 'Pupil', 'Whisk'};
clim = [-5 5];
fh=[];

saveName = 'newCsPlus_image';
fh(end+1) = ensureFigure(saveName, 1);
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
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs+');
set(gcf, 'Position', [304   217   633   485]);

saveName = 'newCsMinus_image';
fh(end+1) = ensureFigure(saveName, 1);
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
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs-');
set(gcf, 'Position', [304   217   633   485]);

saveName = 'alwaysCsPlus_image';
fh(end+1) = ensureFigure(saveName, 1);
cLimFactor = 3;
xData = [min(alwaysCsPlus_trialNumber), max(alwaysCsPlus_trialNumber)];
xlim = [-30 70];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = alwaysCsPlus.(sfield)(sortOrder, :);
    cData = smoothdata(cData, 2, 'movmean', 5, 'omitnan');
    imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:length(sortOrder), [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 length(sortOrder)]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('Always Cs+');
set(gcf, 'Position', [304   217   633   485]);

%% averages
common_odor3 = (-9 <= odor3_trialNumber) & (odor3_trialNumber <= 9);
xlim = [-30 30];
ylim = [-1 2];

% new cs plus
common = sum(~isnan(newCsPlus.licks_cs)) > 3;
savename = 'reversals_newCsPlus';
fh(end + 1) = ensureFigure(savename, 1);
hla = zeros(1,6);
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus.phPeakMean_cs_ch2(goodReversals, common)), nanSEM(newCsPlus.phPeakMean_cs_ch2(goodReversals, common))',...
    'cmap', [237 125 49]/256, 'nan', 'gap'); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus.phPeakMean_cs_ch1(goodReversals, common)), nanSEM(newCsPlus.phPeakMean_cs_ch1(goodReversals, common))',...
    'cmap', [171 55 214]/256, 'nan', 'gap');
hla(2) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus.licks_cs(goodReversals, common)), nanSEM(newCsPlus.licks_cs(goodReversals, common))',...
    'cmap', [0.5 0.5 0.5]/256, 'nan', 'gap');
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch2(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch1(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.licks_cs(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
    '\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\color[rgb]{0.5,0.5,0.5}Odor 3'},...
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');

title('New Cs+');
xlabel('Odor presentations from reversal');
ylabel('Cue response');
formatFigurePoster([5.5 4], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    


% new cs minus
common = sum(~isnan(newCsMinus.licks_cs)) > 3;

savename = 'reversals_newCsMinus';
fh(end + 1) = ensureFigure(savename, 1);
hla = zeros(1,6);
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_ch2(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_ch2(goodReversals, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_ch1(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_ch1(goodReversals, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.licks_cs(goodReversals, common)), nanSEM(newCsMinus.licks_cs(goodReversals, common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch2(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch1(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.licks_cs(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);    
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
            '\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\color[rgb]{0.5,0.5,0.5}Odor 3'},...
            'Location', 'southwest', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
title('New Cs-');
xlabel('Odor presentations from reversal');
ylabel('Cue response');    
formatFigurePoster([5.5 4], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

% always cs plus
common = sum(~isnan(alwaysCsPlus.licks_cs)) > 3;
savename = 'reversals_alwaysCsPlus';
fh(end + 1) = ensureFigure(savename, 1);
hla = zeros(1,6);
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus.phPeakMean_cs_ch2(goodReversals, common)), nanSEM(alwaysCsPlus.phPeakMean_cs_ch2(goodReversals, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus.phPeakMean_cs_ch1(goodReversals, common)), nanSEM(alwaysCsPlus.phPeakMean_cs_ch1(goodReversals, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus.licks_cs(goodReversals, common)), nanSEM(alwaysCsPlus.licks_cs(goodReversals, common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch2(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch1(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.licks_cs(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);    
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
    '\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\color[rgb]{0.5,0.5,0.5}Odor 3'},...
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
title('Always Cs+');
xlabel('Odor presentations from reversal');
ylabel('Cue response');    
formatFigurePoster([5.5 4], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end 

% odor 3
savename = 'reversals_odor3';
ensureFigure(savename, 1);
o3xdata = (0:(size(AR.thirdOdor.licks_cs.before, 2) + size(AR.thirdOdor.licks_cs.after, 2) - 1)) - size(AR.thirdOdor.licks_cs.before, 2) - 1; 
o3licks = [AR.thirdOdor.licks_cs.before AR.thirdOdor.licks_cs.after];
o3_ch1 = [AR.thirdOdor.phPeakMean_cs_ch1.before AR.thirdOdor.phPeakMean_cs_ch1.after];
o3_ch2 = [AR.thirdOdor.phPeakMean_cs_ch2.before AR.thirdOdor.phPeakMean_cs_ch2.after];
hla = zeros(1,3);
[hl, hp] = boundedline(o3xdata, nanmean(odor3.phPeakMean_cs_ch2(goodReversals, :)), nanSEM(odor3.licks_cs(goodReversals, :))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(o3xdata, nanmean(odor3.phPeakMean_cs_ch1(goodReversals, :)), nanSEM(odor3.phPeakPercentile_cs_ch1(goodReversals, :))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(o3xdata, nanmean(odor3.licks_cs(goodReversals, :)), nanSEM(odor3.phPeakMean_cs_ch2(goodReversals, :))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
set(hla, 'LineWidth', 2);    
set(gca, 'XLim', [-10 10]);%, 'YLim', [-1 2]);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks'},...
            'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
title('Odor 3');
xlabel('Odor 3 presentations from reversal');
ylabel('Cue response');
formatFigurePoster([5.5 4], '', 12);
set(gca, 'XLim', [-10 10]);%, 'YLim', [-1 2]);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end


%% reversal averages for whisk and pupil
common_odor3 = (-9 <= odor3_trialNumber) & (odor3_trialNumber <= 9);

% new Cs+
common = sum(~isnan(newCsPlus.licks_cs)) > 3;
xlim = [-30 30];
savename = 'reversals_pupWhisk_newCsPlus';
fh(end + 1) = ensureFigure(savename, 1);
hla = [];
subplot(1,2,1);
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus.pupil_csBaselined(goodReversals & goodPupil, common)), nanSEM(newCsPlus.pupil_csBaselined(goodReversals & goodPupil, common))',...
    'k', 'nan', 'gap'); hold on
hla(end+1) = hl;
hla(end+1) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.pupil_csBaselined(goodReversals & goodPupil, common_odor3)), 'k--');
title('New Cs+');
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'new CS+', 'odor 3'},...        
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
xlabel('Odor presentations from reversal');
ylabel('Pupil cue response (baselined)');

subplot(1,2,2);
hla = [];
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus.whisk_cs(goodReversals & goodWhisk, common)), nanSEM(newCsPlus.whisk_cs(goodReversals & goodWhisk, common))',...
    'k', 'nan', 'gap'); hold on;
hla(end+1) = hl;
hla(end+1) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.whisk_cs(goodReversals & goodWhisk, common_odor3)), 'k--');
title('New Cs+');
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'new CS+', 'odor 3'},...        
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
xlabel('Odor presentations from reversal');
ylabel('Whisk cue response');
formatFigurePoster([10 4], '', 12);

% new Cs-
common = sum(~isnan(newCsMinus.licks_cs)) > 3;
xlim = [-30 30];
savename = 'reversals_pupWhisk_newCsMinus';
fh(end + 1) = ensureFigure(savename, 1);
hla = [];
subplot(1,2,1);
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.pupil_csBaselined(goodReversals & goodPupil, common)), nanSEM(newCsMinus.pupil_csBaselined(goodReversals & goodPupil, common))',...
    'k', 'nan', 'gap'); hold on
hla(end+1) = hl;
hla(end+1) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.pupil_csBaselined(goodReversals & goodPupil, common_odor3)), 'k--');
title('New Cs-');
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'new CS-', 'odor 3'},...        
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
xlabel('Odor presentations from reversal');
ylabel('Pupil cue response (baselined)');

subplot(1,2,2);
hla = [];
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.whisk_cs(goodReversals & goodWhisk, common)), nanSEM(newCsMinus.whisk_cs(goodReversals & goodWhisk, common))',...
    'k', 'nan', 'gap'); hold on;
hla(end+1) = hl;
hla(end+1) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.whisk_cs(goodReversals & goodWhisk, common_odor3)), 'k--');
title('New Cs-');
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'new CS-', 'odor 3'},...        
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
xlabel('Odor presentations from reversal');
ylabel('Whisk cue response');
formatFigurePoster([10 4], '', 12)

% always Cs+
common = sum(~isnan(alwaysCsPlus.licks_cs)) > 3;
xlim = [-30 30];
savename = 'reversals_pupWhisk_alwaysCsPlus';
fh(end + 1) = ensureFigure(savename, 1);
hla = [];
subplot(1,2,1);
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus.pupil_csBaselined(goodReversals & goodPupil, common)), nanSEM(alwaysCsPlus.pupil_csBaselined(goodReversals & goodPupil, common))',...
    'k', 'nan', 'gap'); hold on
hla(end+1) = hl;
hla(end+1) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.pupil_csBaselined(goodReversals & goodPupil, common_odor3)), 'k--');
title('Always Cs+');
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'always CS+', 'odor 3'},...        
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
xlabel('Odor presentations from reversal');
ylabel('Pupil cue response (baselined)');

subplot(1,2,2);
hla = [];
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus.whisk_cs(goodReversals & goodWhisk, common)), nanSEM(alwaysCsPlus.whisk_cs(goodReversals & goodWhisk, common))',...
    'k', 'nan', 'gap'); hold on;
hla(end+1) = hl;
hla(end+1) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.whisk_cs(goodReversals & goodWhisk, common_odor3)), 'k--');
title('Always Cs+');
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'Always CS+', 'odor 3'},...        
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
xlabel('Odor presentations from reversal');
ylabel('Whisk cue response');
formatFigurePoster([10 4], '', 12)

% %% averages, Reweard response, deconvolved
% common_odor3 = (-9 <= odor3_trialNumber) & (odor3_trialNumber <= 9);
% xlim = [-30 30];
% ylim = [-1 2];
% 
% 
% % always cs plus reward
% common = sum(~isnan(alwaysCsPlusReward.licks_cs)) > 3;
% savename = 'reversals_alwaysCsPlusDeconv';
% fh(end + 1) = ensureFigure(savename, 1);
% hla = []'
% [hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlusReward.phPeakMean_us_ch2_deconv(goodReversals, common)), nanSEM(alwaysCsPlusReward.phPeakMean_us_ch2_deconv(goodReversals, common))',...
%     'cmap', [237 125 49]/256); hold on
% hla(end+1) = hl;
% [hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlusReward.phPeakMean_us_ch1_deconv(goodReversals, common)), nanSEM(alwaysCsPlusReward.phPeakMean_us_ch1_deconv(goodReversals, common))',...
%     'cmap', [171 55 214]/256);
% hla(end+1) = hl;
% 
% hla(end+1) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_us_ch2_deconv(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
% hla(end+1) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_us_ch1_deconv(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
% 
% set(hla, 'LineWidth', 2);
% set(gca, 'XLim', xlim);%, 'YLim', ylim);    
% h  = addOrginLines;
% set(h, 'LineWidth', 2);
% % legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
% %     '\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\color[rgb]{0.5,0.5,0.5}Odor 3'},...
% %     'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
% title('Always Cs+');
% xlabel('CS+ and reward presentations from reversal');
% ylabel('Reward response');    
% formatFigurePoster([5.5 4], '', 12);
% 
% if saveOn
%     saveas(gcf, fullfile(savepath, [savename '.fig']));
%     saveas(gcf, fullfile(savepath, [savename '.jpg']));   
%     saveas(gcf, fullfile(savepath, [savename '.epsc']));   
% end 
% 


%% write to pdf

% h = waitbar(0, 'slowly writing pdfs');
% 
% pdfname = fullfile(DB.path, 'pooled', 'reversals_noPunish_pooled.pdf');
% for counter = 1:length(fh)    
%     if counter == 1
%         export_fig(fh(counter),pdfname);  % write to pdf
%     else
%         export_fig(fh(counter),'-append',pdfname);  % write to pdf
%     end
%     waitbar(counter/length(fh));
% end
% close(h);



%% to detect latency to learning, fit weibull function, exponential, also changepoint detection

% weibull doesn't work great due to overfitting and noisy data
baselineTrials = 20;

fitField = 'licks_cs';

fitData = [AR.csMinus.(fitField).before(:, end - baselineTrials + 1:end) AR.csPlus.(fitField).after(:, :)];% - nanmean(AR.csMinus.licks_cs.before(goodReversals, end - baselineTrials + 1:end), 2); % zero/baseline data at start
model = 'a * (1 - exp(-1 * (x/b)^c)) + d'; % weibull function, CDF form
weibull = struct('object', [], 'gof', [], 'output', [], 'toFit', []);
weibull = repmat(weibull, size(fitData, 1), 1);



for counter = 1:size(fitData, 1)
    toFit = fitData(counter, ~isnan(fitData(counter, :)));
    fo = fitoptions('Method', 'NonlinearLeastSquares', 'Robust','On',... 
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
    
%     weibull(counter).xData = xData - baselineTrials;
end

%% fit exponentials to post-reversal data for both acquisition and extinction
% ss = struct('object', zeros(sum(goodReversals), 1), 'gof', zeros(sum(goodReversals), 1), 'output', zeros(sum(goodReversals), 1), 'toFit', zeros(sum(goodReversals), 1));

compFields = {'csPlus', 'csMinus'};
fitFields = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'};
outputFields = {'object', 'gof', 'output', 'toFit', 'a', 'b', 'c'};
expFit = struct();
for counter = 1:length(compFields)
    for counter2 = 1:length(fitFields)
        for counter3 = 1:length(outputFields)
            if ~ismember(outputFields{counter3}, {'a', 'b', 'c'})
                expFit.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = cell(sum(goodReversals), 1);
            else
                expFit.(compFields{counter}).(fitFields{counter2}).(outputFields{counter3}) = NaN(sum(goodReversals), 1);
            end
        end
    end
end

expModel = 'a + b * exp(-x/c)';
% these options should work for newCsPlus and newCsMinus (up or down and
% for z scored and lick rate data (both should fall within interval of -100
% to 100
fo = fitoptions('Method', 'NonlinearLeastSquares',...
    'Upper', [100  100 1000],...
    'Lower', [-100 -100 0],...    
    'StartPoint', [0 1 10]...
    );
for compCounter = 1:length(compFields)
    for fieldCounter = 1:length(fitFields)
        expFitData = AR.(compFields{compCounter}).(fitFields{fieldCounter}).after(goodReversals, :);
        for counter = 1:size(expFitData, 1)        
            toFit = expFitData(counter, ~isnan(expFitData(counter, :)));
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


%% plot cumsum and do changepoint detection

baselineTrials = 20;
cp_licks = bpChangePoints([AR.csMinus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csPlus.licks_cs.after(:, 1:end)], 2, 1000);

%% plot example reversals with weibull fits, changepoints
nShow = 6;
% good revs include [77 3 48 74];,   3 and 74 best examples
% mediocre,  84 and 91 and 42 are mediocre ones
% 4 is good but gradual
toShow = find(goodReversals);

ordering = randperm(sum(goodReversals), nShow);
toShow = toShow(ordering);
toShow2 = 1:sum(goodReversals);
toShow2 = toShow2(ordering);

% nShow = 4;  toShow = [74 4 91 42];
ensureFigure('reversals_pooled_Latency_examples', 1);

for counter = 1:nShow   
    thisRev = toShow(counter);
    thisRev2 = toShow2(counter);
    % weibull
    subplot(nShow, 3, counter*3 - 2);
    plot(weibull(thisRev).toFit, 'g.'); hold on;
    plot(weibull(thisRev).object, 'predfunc'); legend off;
    set(gca, 'XLim', [0 120]);
    % changepoint
    subplot(nShow, 3, counter*3 - 1); hold on;
    scatter(1:length(cp_licks.cumsum{thisRev}), cp_licks.cumsum{thisRev}, 10, cp_licks.logitAll{thisRev}); colormap jet;
    line(repmat(cp_licks.index(thisRev), 1, 2), get(gca, 'YLim')); set(gca, 'XLim', [0 120]);
    textBox(sprintf('Logit=%.2f', cp_licks.logit(thisRev)));
    % exponential
    subplot(nShow, 3, counter*3); hold on;    
    plot(expFit.csPlus.licks_cs.toFit{thisRev2}, 'g.'); hold on;
    plot(expFit.csPlus.licks_cs.object{thisRev2}, 'predfunc'); legend off;
    set(gca, 'XLim', [0 120]);
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
fh(end+1) = ensureFigure(saveName, 1); colormap parula;
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
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs+');
set(gcf, 'Position', [304   217   633   485]);

%% data and images aligned to lick changepoint
fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'licks_cs', 'whisk_cs'};
titles = {'ACh', 'Dop.', 'Licks', 'Whisk'};
clim = [-5 5];
fh=[];

for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,2,fcounter);
    


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

figure;
subplot(2,2,1); plot(x', [y1' y2']); title('not zscored'); xlabel('x'); ylabel('y');
[r,lags] = xcorr(y1,y2, 20, 'unbiased');
subplot(2,2,3); plot(lags, r); ylabel('xcorr'); xlabel('lags');
subplot(2,2,2); plot(x ,zscore(y1), '-', 'LineWidth', 2); hold on; plot(x ,zscore(y2), '--', 'LineWidth', 4);  title('zscored'); xlabel('x'); ylabel('y (zscored)');
[r,lags] = xcorr(y1 - mean(y1) ,y2 - mean(y2), 20, 'unbiased');
subplot(2,2,4); plot(lags, r); ylabel('xcorr'); xlabel('lags');

%% compare strength of licking correlations, ChAT-cre vs DAT-cre
all_Licks = [reshape(newCsPlus.licks_cs(goodReversals, :), numel(newCsPlus.licks_cs(goodReversals, :)), 1); reshape(newCsMinus.licks_cs(goodReversals, :), numel(newCsMinus.licks_cs(goodReversals, :)), 1)]; 
all_ChAT = [reshape(newCsPlus.phPeakMean_cs_ch1(goodReversals, :), numel(newCsPlus.phPeakMean_cs_ch1(goodReversals, :)), 1); reshape(newCsMinus.phPeakMean_cs_ch1(goodReversals, :), numel(newCsMinus.phPeakMean_cs_ch1(goodReversals, :)), 1)]; 
all_DAT = [reshape(newCsPlus.phPeakMean_cs_ch2(goodReversals, :), numel(newCsPlus.phPeakMean_cs_ch2(goodReversals, :)), 1); reshape(newCsMinus.phPeakMean_cs_ch2(goodReversals, :), numel(newCsMinus.phPeakMean_cs_ch2(goodReversals, :)), 1)]; 
keep = isfinite(all_Licks) & isfinite(all_ChAT) & isfinite(all_DAT);
all_Licks = all_Licks(keep);
all_ChAT = all_ChAT(keep);
all_DAT = all_DAT(keep);

ensureFigure('test_corr', 1); 
subplot(1,2,1); scatter(all_Licks, all_ChAT, 8, '.'); ylabel('cue ChAT'); textBox(sprintf('R=%.2f', corr(all_Licks, all_ChAT)));
subplot(1,2,2); scatter(all_Licks, all_DAT, 8, '.'); ylabel('cue DAT'); textBox(sprintf('R=%.2f', corr(all_Licks, all_DAT)));


