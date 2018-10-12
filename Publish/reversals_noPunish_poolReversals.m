
DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, 'pooled', filesep);
smoothWindow = 1;
saveOn = 1;
%% create structure containing all reversals
AR = struct('csPlus', [], 'csMinus', [], 'csPlusReward', [], 'thirdOdor', []); % AR = all reversals
for si = 1:length(DB.animals)    
    animal = DB.animals{si};
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
trialWindow = [-30 30];%?

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
rocThresh = 0.6;
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

%% filter reversals according to quality

goodReversals = ...
    ~isnan(trialsToCriterion);% &...
%     auROC.phPeakMean_cs_ch1.before > 0.2 &...
%     auROC.phPeakMean_cs_ch2.before > 0.2;
sortVariable = trialsToCriterion;
[~, sortOrder] = sort(sortVariable);
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

saveName = 'newCsPlus_image';
fh(end+1) = ensureFigure(saveName, 1);
cLimFactor = 3;
xlim = [min(newCsPlus_trialNumber), max(newCsPlus_trialNumber)];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = newCsPlus.(sfield)(sortOrder, :);
    imagesc('XData', xlim, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 nReversals]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs+');
set(gcf, 'Position', [304   217   633   485]);

saveName = 'newCsMinus_image';
fh(end+1) = ensureFigure(saveName, 1);
cLimFactor = 3;
xlim = [min(newCsMinus_trialNumber), max(newCsMinus_trialNumber)];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = newCsMinus.(sfield)(sortOrder, :);
    imagesc('XData', xlim, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 nReversals]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs-');
set(gcf, 'Position', [304   217   633   485]);

saveName = 'alwaysCsPlus_image';
fh(end+1) = ensureFigure(saveName, 1);
cLimFactor = 3;
xlim = [min(alwaysCsPlus_trialNumber), max(alwaysCsPlus_trialNumber)];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = alwaysCsPlus.(sfield)(sortOrder, :);
    imagesc('XData', xlim, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 nReversals]);
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

h = waitbar(0, 'slowly writing pdfs');

pdfname = fullfile(DB.path, 'pooled', 'reversals_noPunish_pooled.pdf');
for counter = 1:length(fh)    
    if counter == 1
        export_fig(fh(counter),pdfname);  % write to pdf
    else
        export_fig(fh(counter),'-append',pdfname);  % write to pdf
    end
    waitbar(counter/length(fh));
end
close(h);




