
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
% goodReversals = ...
%     ~isnan(trialsToCriterion) &...
%     auROC.phPeakMean_cs_ch1.acq > 0 &...
%     auROC.phPeakMean_cs_ch2.acq > 0;

goodReversals = ...
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
    'phPeakMean_us_ch1_deconv', 'phPeakMean_us_ch2_deconv', 'phPeakPercentile_us_ch1_deconv', 'phPeakPercentile_us_ch2_deconv', 'phBaseline_ch1', 'phBaseline_ch2', 'phPeakMean_us_ch1', 'phPeakMean_us_ch2'};

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

newCsPlus_trialNumber = (1:size(newCsPlus.licks_cs, 2)) - size(AR.csMinus.licks_cs.before, 2); newCsPlus.trialNumber = newCsPlus_trialNumber;
newCsPlus_firstRevTrial = size(AR.csMinus.licks_cs.before, 2) + 1; newCsPlus.firstRevTrial = newCsPlus_firstRevTrial;
newCsMinus_trialNumber = (1:size(newCsMinus.licks_cs, 2)) - size(AR.csPlus.licks_cs.before, 2);  newCsMinus.trialNumber = newCsMinus_trialNumber;
newCsMinus_firstRevTrial = size(AR.csPlus.licks_cs.before, 2) + 1; newCsMinus.firstRevTrial = newCsMinus_firstRevTrial;
alwaysCsPlus_trialNumber = (1:size(alwaysCsPlus.licks_cs, 2)) - size(AR.csPlus.licks_cs.before, 2); alwaysCsPlus.trialNumber = alwaysCsPlus_trialNumber;
alwaysCsPlus_firstRevTrial = size(AR.csPlus.licks_cs.before, 2) + 1; alwaysCsPlus.firstRevTrial = alwaysCsPlus_firstRevTrial;
alwaysCsPlusReward_trialNumber = (1:size(alwaysCsPlusReward.licks_cs, 2)) - size(AR.csPlusReward.licks_cs.before, 2); alwaysCsPlusReward.trialNumber = alwaysCsPlusReward_trialNumber;
alwaysCsPlusReward_firstRevTrial = size(AR.csPlusReward.licks_cs.before, 2) + 1; alwaysCsPlusReward.firstRevTrial = alwaysCsPlusReward_firstRevTrial; 
odor3_trialNumber = (1:size(odor3.licks_cs, 2)) - size(AR.thirdOdor.licks_cs.before, 2); 

oldCsPlus_trialNumber = -(size(AR.csPlus.licks_cs.before, 2) - 1) : 0;
oldCsMinus_trialNumber = -(size(AR.csMinus.licks_cs.before, 2) - 1) : 0;

%% take advantage of paired recordings, subtract dopamine cue response from Ach. cue response, etc.
newCsPlus.phPeakMean_cs_AchMinusDop = newCsPlus.phPeakMean_cs_ch1 - newCsPlus.phPeakMean_cs_ch2;
newCsMinus.phPeakMean_cs_AchMinusDop = newCsMinus.phPeakMean_cs_ch1 - newCsMinus.phPeakMean_cs_ch2;
alwaysCsPlus.phPeakMean_cs_AchMinusDop = alwaysCsPlus.phPeakMean_cs_ch1 - alwaysCsPlus.phPeakMean_cs_ch2;
odor3.phPeakMean_cs_AchMinusDop = odor3.phPeakMean_cs_ch1 - odor3.phPeakMean_cs_ch2;

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

%% scatter plots (alternative to averages)
common_odor3 = (-9 <= odor3_trialNumber) & (odor3_trialNumber <= 9);

ylim = [-1 2];
markerSize = 5;
marker = '.';

% new cs plus
common = sum(~isnan(newCsPlus.licks_cs)) > 3;
savename = 'reversals_newCsPlus_scatter';
fh(end + 1) = ensureFigure(savename, 1);

subplot(2,2,1);
dataX = repmat(newCsPlus_trialNumber(common), sum(goodReversals), 1); dataX = dataX(:);
dataY = newCsPlus.phPeakMean_cs_ch2(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, [237 125 49]/256, marker); hold on
h  = addOrginLines; 
ylabel('Dop.');
title('New Cs+');

subplot(2,2,2);
dataY = newCsPlus.phPeakMean_cs_ch1(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, [171 55 214]/256, marker); hold on
h  = addOrginLines; 
xlabel('Odor presentations from rev.'); 
ylabel('ACh.');

subplot(2,2,3);
dataY = newCsPlus.licks_cs(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, [0.5 0.5 0.5]/256, marker); hold on
h  = addOrginLines; 
xlabel('Odor presentations from rev.'); 
ylabel('Licks');
formatFigurePublish('size', [3 2]);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

% new cs minus
common = sum(~isnan(newCsMinus.licks_cs)) > 3;
savename = 'reversals_newCsMinus_scatter';
fh(end + 1) = ensureFigure(savename, 1);

subplot(2,2,1);
dataX = repmat(newCsMinus_trialNumber(common), sum(goodReversals), 1); dataX = dataX(:);
dataY = newCsMinus.phPeakMean_cs_ch2(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, [237 125 49]/256, marker); hold on
h  = addOrginLines; 
ylabel('Dop.');
title('New Cs-');

subplot(2,2,2);
dataY = newCsMinus.phPeakMean_cs_ch1(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, [171 55 214]/256, marker); hold on
h  = addOrginLines; 
xlabel('Odor presentations from rev.'); 
ylabel('ACh.');

subplot(2,2,3);
dataY = newCsMinus.licks_cs(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, [0.5 0.5 0.5]/256, marker); hold on
h  = addOrginLines; 
xlabel('Odor presentations from rev.'); 
ylabel('Licks');
formatFigurePublish('size', [3 2]);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

% always cs plus
common = sum(~isnan(alwaysCsPlus.licks_cs)) > 3;
savename = 'reversals_alwaysCsPlus_scatter';
fh(end + 1) = ensureFigure(savename, 1);

subplot(2,2,1);
dataX = repmat(alwaysCsPlus_trialNumber(common), sum(goodReversals), 1); dataX = dataX(:);
dataY = alwaysCsPlus.phPeakMean_cs_ch2(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, [237 125 49]/256, marker); hold on
h  = addOrginLines; 
ylabel('Dop.');
title('Always Cs+');

subplot(2,2,2);
dataY = alwaysCsPlus.phPeakMean_cs_ch1(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, [171 55 214]/256, marker); hold on
h  = addOrginLines; 
xlabel('Odor presentations from rev.'); 
ylabel('ACh.');

subplot(2,2,3);
dataY = alwaysCsPlus.licks_cs(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, [0.5 0.5 0.5]/256, marker); hold on
h  = addOrginLines; 
xlabel('Odor presentations from rev.'); 
ylabel('Licks');
formatFigurePublish('size', [3 2]);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

%% test subtraction approach to take advantage of paired recordings

savename = 'ACh_minus_Dop_avgs';
ensureFigure(savename, 1);
trialRange = [-30 30];
subplot(2,2,1);
title('new Cs+');
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common)), nanSEM(newCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common))',...
    'cmap', [0 1 0], 'nan', 'gap'); hold on
set(gca, 'XLim', trialRange);%, 'YLim', [-1 2]);
t = textBox('\color[rgb]{0,0.5,0}ACh. - Dop.'); 
set(t, 'Interpreter', 'tex', 'FontSize', 8, 'FontWeight', 'bold');
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');
addOrginLines;
% t = textBox('\color{green} ACh. minus Dop.'); set(t, 'Interpreter', 'tex', 'FontSize', 8);

subplot(2,2,2);
title('new Cs-');
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_AchMinusDop(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_AchMinusDop(goodReversals, common))',...
    'cmap', [0 1 0], 'nan', 'gap'); hold on
set(gca, 'XLim', trialRange);%, 'YLim', [-1 2]);
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');
addOrginLines;
% t = textBox('\color{green} ACh. minus Dop.'); set(t, 'Interpreter', 'tex', 'FontSize', 8);
subplot(2,2,3);
title('always Cs+');
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common)), nanSEM(alwaysCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common))',...
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


%% scatter version of subtraction approach to take advantage of paired recordings

xlim = [-30 30];
ylim = [-1 2];
markerSize = 5;
marker = '.';
rc = [0 0.75 0];

savename = 'ACh_minus_Dop_scatter';
ensureFigure(savename, 1);
subplot(2,2,1); hold on;
title('new Cs+');
dataX = repmat(newCsPlus_trialNumber(common), sum(goodReversals), 1); dataX = dataX(:);
dataY = newCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, rc, marker); hold on
t = textBox('\color[rgb]{0,0.5,0}ACh. - Dop.'); 
set(t, 'Interpreter', 'tex', 'FontSize', 8, 'FontWeight', 'bold');


ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');
addOrginLines;
% t = textBox('\color{green} ACh. minus Dop.'); set(t, 'Interpreter', 'tex', 'FontSize', 8);

subplot(2,2,2); hold on;
title('new Cs-');
dataX = repmat(newCsMinus_trialNumber(common), sum(goodReversals), 1); dataX = dataX(:);
dataY = newCsMinus.phPeakMean_cs_AchMinusDop(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, rc, marker); hold on


ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');
addOrginLines;
% t = textBox('\color{green} ACh. minus Dop.'); set(t, 'Interpreter', 'tex', 'FontSize', 8);

subplot(2,2,3); hold on;
title('always Cs+');
dataX = repmat(alwaysCsPlus_trialNumber(common), sum(goodReversals), 1); dataX = dataX(:);
dataY = alwaysCsPlus.phPeakMean_cs_AchMinusDop(goodReversals, common); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, rc, marker); hold on
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)'); xlabel('trials from rev.');
addOrginLines;
% t = textBox('\color{green} ACh. minus Dop.'); set(t, 'Interpreter', 'tex', 'FontSize', 8);

subplot(2,2,4); hold on;
title('odor 3');
dataX = repmat(odor3_trialNumber(common_odor3), sum(goodReversals), 1); dataX = dataX(:);
dataY = odor3.phPeakMean_cs_AchMinusDop(goodReversals, common_odor3); dataY = dataY(:);
hl = scatter(dataX, dataY, markerSize, rc, marker); hold on
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)'); xlabel('trials from rev.');
addOrginLines;


% subplot(2,2,4);
% [ax, t] = textAxes(gca, {'\color{green}Green=', '\color{green}ACh. - Dop.'}); 
% set(t, 'Interpreter', 'tex', 'FontSize', 8, 'FontWeight', 'bold');
formatFigurePublish('size', [3 2]);
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



%% to detect latency to learning, fit weibull function, exponential functions, also changepoint detection
baselineTrials = 20;
fitField = 'phPeakMean_cs_ch1';

%%
% %% fit weibull function
% 
% % weibull doesn't work great due to overfitting and noisy data
% 
% 
% fitData = [AR.csMinus.(fitField).before(:, end - baselineTrials + 1:end) AR.csPlus.(fitField).after(:, :)];% - nanmean(AR.csMinus.licks_cs.before(goodReversals, end - baselineTrials + 1:end), 2); % zero/baseline data at start
% model = 'a * (1 - exp(-1 * (x/b)^c)) + d'; % weibull function, CDF form
% weibull = struct('object', [], 'gof', [], 'output', [], 'toFit', []);
% weibull = repmat(weibull, size(fitData, 1), 1);
% 
% 
% 
% for counter = 1:size(fitData, 1)
%     toFit = fitData(counter, ~isnan(fitData(counter, :)));
%     fo = fitoptions('Method', 'NonlinearLeastSquares',... 
%         'Upper', [Inf  Inf Inf Inf],...  % 20 (3rd upper)
%         'Lower', [0 0 0 -Inf],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
%         'StartPoint', [mean(toFit) baselineTrials baselineTrials min(toFit)]...
%         );
%     
%     ft = fittype(model, 'options', fo);
% %     xData = (0:length(toFit) - 1)';
%     [fitobject, gof, output] = fit((0:length(toFit) - 1)', toFit', ft, fo);
%     weibull(counter).object = fitobject;
%     weibull(counter).gof = gof;
%     weibull(counter).output = output;
%     weibull(counter).toFit = toFit;
%     
% %     weibull(counter).xData = xData - baselineTrials;
% end


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



%% changepoint detection for acquisition (acq) and extinction (ext)
% compFields = {'csPlus', 'csMinus'};
% fitFields = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'};


% cp_licks = bpChangePoints([AR.csMinus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csPlus.licks_cs.after(:, 1:end)], 2, 1000);
cp.csPlus.licks_cs = bpChangePoints([AR.csMinus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csPlus.licks_cs.after(:, 1:end)], 2, 1000);
cp.csPlus.phPeakMean_cs_ch1 = bpChangePoints([AR.csMinus.phPeakMean_cs_ch1.before(:, end - baselineTrials + 1:end) AR.csPlus.phPeakMean_cs_ch1.after(:, 1:end)], 2, 1000);
cp.csPlus.phPeakMean_cs_ch2 = bpChangePoints([AR.csMinus.phPeakMean_cs_ch2.before(:, end - baselineTrials + 1:end) AR.csPlus.phPeakMean_cs_ch2.after(:, 1:end)], 2, 1000);
cp.csMinus.licks_cs = bpChangePoints([AR.csPlus.licks_cs.before(:, end - baselineTrials + 1:end) AR.csMinus.licks_cs.after(:, 1:end)], 2, 1000);
cp.csMinus.phPeakMean_cs_ch1 = bpChangePoints([AR.csPlus.phPeakMean_cs_ch1.before(:, end - baselineTrials + 1:end) AR.csMinus.phPeakMean_cs_ch1.after(:, 1:end)], 2, 1000);
cp.csMinus.phPeakMean_cs_ch2 = bpChangePoints([AR.csPlus.phPeakMean_cs_ch2.before(:, end - baselineTrials + 1:end) AR.csMinus.phPeakMean_cs_ch2.after(:, 1:end)], 2, 1000);

%% make a bar graph or box plot of changepoints
    
savename = 'changepoints_all';
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
    ydata = cp.(fields{counter, 1}).(fields{counter, 2}).index(goodReversals) - baselineTrials;   
    kstest(ydata)
    all_cps = [all_cps ydata];
    scatter(repmat(counter, numel(ydata), 1) + (rand(numel(ydata), 1) - 0.5)/3, ydata, markerSize, fields{counter, 3}, '.');
    errorbar(counter, mean(ydata), std(ydata)/sqrt(numel(ydata)), 'Color', fields{counter, 3}, 'LineWidth', 2)
end

set(gca, 'XLim', [0.5 6.5], 'Ylim', [-20 40], 'XTick', 1:6, 'XTickLabel', fields(:,4)', 'FontSize', 12); ylabel('trials from rev.', 'FontSize', 12);


formatFigurePublish('size', [3.5 2], 'fontSize', 12);
if saveOn 
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    export_fig(fullfile(savepath, savename), '-eps');
end

%% despite paired measurements, weak correlations between detected changepoints
savename = 'ChangePoint_correlations_scatter';
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

%% (DOESN'T LOOK GOOD) make a bar graph of changepoints, this time connect the lines between the paired conditions
    
% saveName = 'changepoints_all_conected';
% ensureFigure(savename, 1); axes('FontSize', 12); hold on;
% fields = {'licks_acq', 'chat_acq', 'dat_acq', 'licks_ext', 'chat_ext', 'dat_ext'};
% colors = {mycolors('licks') mycolors('chat') mycolors('dat') mycolors('licks') mycolors('chat') mycolors('dat')};
% 
% all_cps = [];
% for counter = 1:length(fields)
%     ydata = cp.(fields{counter}).index(goodReversals) - baselineTrials;
%     all_cps = [all_cps ydata];
% end
% 
% plot([1:3 NaN 4:6]', [all_cps(:,1:3) NaN(size(all_cps, 1), 1) all_cps(:,4:6)]', '-', 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]);
% 
% errorbar((1:3), mean(all_cps(:,1:3)), std(all_cps(:,1:3))/sqrt(size(all_cps, 1)), 'Color', colors{counter}, 'LineWidth', 2)
% errorbar((4:6), mean(all_cps(:,4:6)), std(all_cps(:,4:6))/sqrt(size(all_cps, 1)), 'Color', colors{counter}, 'LineWidth', 2)
% 
% set(gca, 'XLim', [0.5 6.5], 'Ylim', [-20 40], 'XTick', 1:6, 'XTickLabel', {'licks', 'Ach.', 'Dop.', 'licks', 'Ach.', 'Dop.'}, 'FontSize', 12); ylabel('trials from rev.', 'FontSize', 12);
% 
% 
% formatFigurePublish('size', [3.5 2], 'fontSize', 12);
% if saveOn 
%     export_fig(fullfile(savepath, saveName), '-eps');
% end

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
% xData = [min(newCsPlus_trialNumber), max(newCsPlus_trialNumber)];
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

