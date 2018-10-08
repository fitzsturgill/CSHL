
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
comp = {'licks_cs', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'}; % comparisons
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

%% filter reversals according to quality

goodReversals = ...
    ~isnan(trialsToCriterion);% &...
%     auROC.phPeakMean_cs_ch1.before > 0.2 &...
%     auROC.phPeakMean_cs_ch2.before > 0.2;
sortVariable = trialsToCriterion;
[~, sortOrder] = sort(sortVariable);
%

%% compile data
fieldsToCompile = {...
    'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'csLicksROC', 'licks_cs', 'pupil_cs', 'whisk_cs', 'wheel_baseline'};

newCsPlus = struct();
newCsMinus = struct();
alwaysCsPlus = struct();
odor3 = struct();

newCsPlus.trialNumber = (1:size(newCsPlus_ch1, 2)) - size(AR.csMinus.phPeakMean_cs_ch1.before, 2);
newCsPlus.firstRevTrial = size(AR.csMinus.phPeakMean_cs_ch1.before, 2) + 1;
newCsMinus.trialNumber = (1:size(newCsMinus_ch1, 2)) - size(AR.csPlus.phPeakMean_cs_ch1.before, 2);
newCsMinus.firstRevTrial = size(AR.csPlus.phPeakMean_cs_ch1.before, 2) + 1;
alwaysCsPlus.trialNumber = (1:size(alwaysCsPlus_ch1, 2)) - size(AR.csPlus.phPeakMean_cs_ch1.before, 2);
alwaysCsPlus.firstRevTrial = size(AR.csPlus.phPeakMean_cs_ch1.before, 2) + 1;
odor3.trialNumber = (1:size(odor3_ch1, 2)) - size(AR.thirdOdor.phPeakMean_cs_ch1.before, 2);

for counter = 1:length(fieldToCompile)
    field = fieldsToCompile{counter};
    newCsPlus.(field) = smoothdata([AR.csMinus.(field).before AR.csPlus.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    newCsMinus.(field) = smoothdata([AR.csPlus.(field).before AR.csMinus.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    alwaysCsPlus.(field) = smoothdata([AR.csPlus.(field).before AR.csPlus.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
    odor3.(field) = smoothdata([AR.thirdOdor.(field).before AR.thirdOdor.(field).after], 2, 'movmean', smoothWindow, 'omitnan');
end

oldCsPlus_trialNumber = -(size(AR.csPlus.phPeakMean_cs_ch1.before, 2) - 1) : 0;
oldCsMinus_trialNumber = -(size(AR.csMinus.phPeakMean_cs_ch1.before, 2) - 1) : 0;



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


orderings = {'newCsPlus', 'newCsMinus', 'alwaysCsPlus'};
% fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'csLicksROC', 'licks_cs', 'pupil_cs', 'whisk_cs', 'wheel_baseline'};
clim = [-5 5];
for oix = 1:length(orderings)
    revOrder = orderings{oix};
    saveName = [revOrder '_image'];
    ensureFigure(saveName, 1);
    
    subplot(2,3,1);    
    xlim = [min(newCsPlus_trialNumber), max(newCsPlus_trialNumber)];
    imagesc('XData', xlim, 'CData', newCsPlus_ch1_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('ACh');  set(gca, 'CLim', clim)
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 nReversals]);
    ylabel('Reversal #');
    subplot(2,2,2);
    imagesc('XData', xlim, 'CData', newCsPlus_ch2_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Dop');  set(gca, 'CLim', clim)
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
    set(gca, 'YLim', [1 nReversals]);
    subplot(2,2,3);
    imagesc('XData', xlim, 'CData', newCsPlus_licks_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on;title('Licks');  set(gca, 'CLim', clim)
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
    set(gca, 'YLim', [1 nReversals]);
    xlabel('new Cs+ trials from reversal');
    ylabel('Reversal #');
    subplot(2,2,4);
    imagesc('XData', xlim, 'CData', newCsPlus_roc(sortOrder, :)); set(gca, 'XLim', xlim); hold on;title('Lick auROC');  %set(gca, 'CLim', clim)
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
    set(gca, 'YLim', [1 nReversals]);
    xlabel('new Cs+ trials from reversal');
    set(gcf, 'Position', [304   217   633   485]);
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
        disp('figure saved');
    end


% newCsMinus
saveName = 'newCsMinus_image';
ensureFigure(saveName, 1);
subplot(2,2,1);
xlim = [min(newCsMinus_trialNumber), max(newCsMinus_trialNumber)];
imagesc('XData', xlim, 'CData', newCsMinus_ch1_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Ach');  set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
ylabel('Reversal #');
subplot(2,2,2);
imagesc('XData', xlim, 'CData', newCsMinus_ch2_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Dop');  set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
subplot(2,2,3);
imagesc('XData', xlim, 'CData', newCsMinus_licks_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Licks');  set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
xlabel('new Cs- trials from reversal');
ylabel('Reversal #');
subplot(2,2,4);
imagesc('XData', xlim, 'CData', newCsMinus_roc(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Lick auROC');  %set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
xlabel('new Cs- trials from reversal');
set(gcf, 'Position', [304   217   633   485]);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end


% alwaysCsPlus 
saveName = 'alwaysCsPlus_image';
ensureFigure(saveName, 1);
subplot(2,2,1);
xlim = [min(alwaysCsPlus_trialNumber), max(alwaysCsPlus_trialNumber)];
imagesc('XData', xlim, 'CData', alwaysCsPlus_ch1_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Ach');  set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
ylabel('Reversal #');
subplot(2,2,2);
imagesc('XData', xlim, 'CData', alwaysCsPlus_ch2_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Dop');  set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
subplot(2,2,3);
imagesc('XData', xlim, 'CData', alwaysCsPlus_licks_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Licks');  set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
ylabel('Reversal #');
xlabel('always Cs+ trials from reversal');
subplot(2,2,4);
imagesc('XData', xlim, 'CData', alwaysCsPlus_roc(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Lick auROC');  %set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
xlabel('always Cs+ trials from reversal');
set(gcf, 'Position', [304   217   633   485]);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end

%  images, odor 3
xlim = [min(odor3_trialNumber), max(odor3_trialNumber)];
saveName = 'odor3_image';
ensureFigure(saveName, 1);
subplot(2,2,1); title('ACh');
imagesc('XData', xlim, 'CData', odor3_ch1_norm(sortOrder, :)); set(gca, 'XLim', xlim); set(gca, 'CLim', clim); hold on;
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
ylabel('Reversal #');
subplot(2,2,2); title('Dop');
imagesc('XData', xlim, 'CData', odor3_ch2_norm(sortOrder, :)); set(gca, 'XLim', xlim); set(gca, 'CLim', clim); hold on;
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
xlabel('Odor 3 trials from reversal');
subplot(2,2,3); title('Licks');
imagesc('XData', xlim, 'CData', odor3_licks_norm(sortOrder, :)); set(gca, 'XLim', xlim); set(gca, 'CLim', clim); hold on;
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
set(gca, 'YLim', [1 nReversals]);
xlabel('Odor 3 trials from reversal');
ylabel('Reversal #');
set(gcf, 'Position', [304   217   633   485]);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end

%% averages
common = sum(~isnan(newCsPlus_ch1_norm)) > 3;
common_odor3 = (-9 <= odor3_trialNumber) & (odor3_trialNumber <= 9);

xlim = [-30 30];
ylim = [-1 2];

savename = 'reversals_newCsPlus';
ensureFigure(savename, 1);
hla = zeros(1,6);
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_ch2_norm(goodReversals, common)), nanSEM(newCsPlus_ch2_norm(goodReversals, common))',...
    'cmap', [237 125 49]/256, 'nan', 'gap'); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_ch1_norm(goodReversals, common)), nanSEM(newCsPlus_ch1_norm(goodReversals, common))',...
    'cmap', [171 55 214]/256, 'nan', 'gap');
hla(2) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_licks_norm(goodReversals, common)), nanSEM(newCsPlus_licks_norm(goodReversals, common))',...
    'cmap', [0.5 0.5 0.5]/256, 'nan', 'gap');
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch1_norm(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch2_norm(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_licks_norm(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
    set(gca, 'XLim', xlim);%, 'YLim', ylim);
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
        '\bf\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\bf\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\bf\color[rgb]{0.5,0.5,0.5}Odor 3'},...
        'Location', 'southeast', 'FontSize', 14, 'Interpreter', 'tex', 'Box', 'off');


    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');

    formatFigurePoster([5.5 4], '', 16);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    
%
% newCsMinus - normalize by f% of pre reversal values, smooth, find first and last common points across reversals






% common = find(mean(newCsMinus_ch1_norm));% applying mean will reveal NaNs
newCsMinus_common = sum(~isnan(newCsMinus_ch1_norm)) > 3;% applying mean will reveal NaNs
newCsMinus_common = sum(~isnan(newCsPlus_ch1_norm))  > 3;

savename = 'reversals_newCsMinus';
ensureFigure(savename, 1);
hla = zeros(1,6);
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmean(newCsMinus_ch2_norm(goodReversals, newCsMinus_common)), nanSEM(newCsMinus_ch2_norm(goodReversals, newCsMinus_common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmean(newCsMinus_ch1_norm(goodReversals, newCsMinus_common)), nanSEM(newCsMinus_ch1_norm(goodReversals, newCsMinus_common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmean(newCsMinus_licks_norm(goodReversals, newCsMinus_common)), nanSEM(newCsMinus_licks_norm(goodReversals, newCsMinus_common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch1_norm(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch2_norm(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_licks_norm(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
    set(gca, 'XLim', xlim);%, 'YLim', ylim);    
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
                '\bf\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\bf\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\bf\color[rgb]{0.5,0.5,0.5}Odor 3'},...
                'Location', 'southwest', 'FontSize', 14, 'Interpreter', 'tex', 'Box', 'off');

    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');    
    formatFigurePoster([5.5 4], '', 16);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    



%
savename = 'reversals_odor3';
ensureFigure(savename, 1);
o3xdata = (0:(size(AR.thirdOdor.licks_cs.before, 2) + size(AR.thirdOdor.licks_cs.after, 2) - 1)) - size(AR.thirdOdor.licks_cs.before, 2) - 1; 
o3licks = [AR.thirdOdor.licks_cs.before AR.thirdOdor.licks_cs.after];
o3_ch1 = [AR.thirdOdor.phPeakMean_cs_ch1.before AR.thirdOdor.phPeakMean_cs_ch1.after];
o3_ch2 = [AR.thirdOdor.phPeakMean_cs_ch2.before AR.thirdOdor.phPeakMean_cs_ch2.after];
hla = zeros(1,3);
[hl, hp] = boundedline(o3xdata, nanmean(odor3_licks_norm(goodReversals, :)), nanSEM(odor3_licks_norm(goodReversals, :))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(o3xdata, nanmean(odor3_ch1_norm(goodReversals, :)), nanSEM(odor3_ch1_norm(goodReversals, :))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(o3xdata, nanmean(odor3_ch2_norm(goodReversals, :)), nanSEM(odor3_ch2_norm(goodReversals, :))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;

set(hla, 'LineWidth', 2);    
    set(gca, 'XLim', [-10 10]);%, 'YLim', [-1 2]);
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks'},...
                'Location', 'southeast', 'FontSize', 14, 'Interpreter', 'tex', 'Box', 'off');

    xlabel('Odor 3 presentations from reversal');
    ylabel('Cue response (norm.)');

    formatFigurePoster([5.5 4], '', 16);
    set(gca, 'XLim', [-10 10]);%, 'YLim', [-1 2]);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    



%

% common = find(mean(alwaysCsPlus_ch1_norm));% applying mean will reveal NaNs
common = sum(~isnan(alwaysCsPlus_ch1_norm)) > 3;% applying mean will reveal NaNs

savename = 'reversals_alwaysCsPlus';
ensureFigure(savename, 1);
hla = zeros(1,6);
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus_ch2_norm(goodReversals, common)), nanSEM(alwaysCsPlus_ch2_norm(goodReversals, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus_ch1_norm(goodReversals, common)), nanSEM(alwaysCsPlus_ch1_norm(goodReversals, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus_licks_norm(goodReversals, common)), nanSEM(alwaysCsPlus_licks_norm(goodReversals, common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch1_norm(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch2_norm(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_licks_norm(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
    set(gca, 'XLim', xlim);%, 'YLim', ylim);    
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
        '\bf\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\bf\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\bf\color[rgb]{0.5,0.5,0.5}Odor 3'},...
        'Location', 'southeast', 'FontSize', 14, 'Interpreter', 'tex', 'Box', 'off');

    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');    
    formatFigurePoster([5.5 4], '', 16);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    


% %% find good example reversals for Ach, Dop, and licking
% sb = repmat(ceil(sqrt(nReversals)), 1, 2);
% savename1 = 'newCsPlus_all_tiled';
% h1 = ensureFigure(savename1, 1);
% savename2 = 'newCsMinus_all_tiled';
% h2 = ensureFigure(savename2, 1);
% for counter = 1:nReversals
%     a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); hold(a1, 'on');
%     plot(newCsPlus_trialNumber, newCsPlus_ch1_norm(counter, :), 'g', 'Parent', a1); 
%     plot(newCsPlus_trialNumber, newCsPlus_ch2_norm(counter, :), 'r', 'Parent', a1);
%     plot(newCsPlus_trialNumber, newCsPlus_licks_norm(counter, :), 'k', 'Parent', a1);    
%     a2 = subplot(sb(1),sb(2),counter, 'Parent', h2); hold(a2, 'on');
%     plot(newCsMinus_trialNumber, newCsMinus_ch1_norm(counter, :), 'g', 'Parent', a2); 
%     plot(newCsMinus_trialNumber, newCsMinus_ch2_norm(counter, :), 'r', 'Parent', a2);
%     plot(newCsMinus_trialNumber, newCsMinus_licks_norm(counter, :), 'k', 'Parent', a2);        
% end
% if saveOn
%     saveas(h1, fullfile(savepath, [savename1 '.fig']));
%     saveas(h1, fullfile(savepath, [savename1 '.jpg']));   
%     saveas(h1, fullfile(savepath, [savename1 '.epsc']));   
%     saveas(h2, fullfile(savepath, [savename2 '.fig']));
%     saveas(h2, fullfile(savepath, [savename2 '.jpg']));   
%     saveas(h2, fullfile(savepath, [savename2 '.epsc']));   
% end    

% %% plot good example reversals for both conditions  (CS- -> CS+,  CS+ -> CS-) 
% ex1 = 11;
% ex2 = 2; %3;
% acolor = [0.6680,0.2148,0.8359];
% dcolor = [0.9258,0.4883,0.1914];
% lcolor = [0.3 0.3 0.3];
% lwidth = 2;
% savename = 'newCsPlus_example';
% smoothFactor = 5;
% smoothMethod = 'movmean';
% ensureFigure(savename, 1); axes; hold on;
% plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch1_norm(ex1, :), smoothMethod, smoothFactor), 'Color', acolor, 'LineWidth', lwidth);
% plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch2_norm(ex1, :), smoothMethod, smoothFactor), 'Color', dcolor, 'LineWidth', lwidth);
% plot(newCsPlus_trialNumber, smoothdata(newCsPlus_licks_norm(ex1, :), smoothMethod, smoothFactor), 'Color', lcolor, 'LineWidth', lwidth);
% set(gca, 'XLim', [-40 40], 'XTick', [], 'YLim', [-1 1.5]);
% addOrginLines;
% formatFigurePoster([5.5 2], '', 20);
% if saveOn
%     saveas(gcf, fullfile(savepath, [savename '.fig']));
%     saveas(gcf, fullfile(savepath, [savename '.jpg']));   
%     saveas(gcf, fullfile(savepath, [savename '.epsc']));   
% end    
% 
% savename = 'newCsMinus_example';
% ensureFigure(savename, 1); axes; hold on;
% plot(newCsMinus_trialNumber, smoothdata(newCsMinus_ch1_norm(ex2, :), smoothMethod, smoothFactor), 'Color', acolor, 'LineWidth', lwidth);
% plot(newCsMinus_trialNumber, smoothdata(newCsMinus_ch2_norm(ex2, :), smoothMethod, smoothFactor), 'Color', dcolor, 'LineWidth', lwidth);
% plot(newCsMinus_trialNumber, smoothdata(newCsMinus_licks_norm(ex2, :), smoothMethod, smoothFactor), 'Color', lcolor, 'LineWidth', lwidth);
% set(gca, 'XLim', [-40 40], 'XTick', [], 'YLim', [-1 1.5]);
% addOrginLines;
% formatFigurePoster([5.5 2], '', 20);
% if saveOn
%     saveas(gcf, fullfile(savepath, [savename '.fig']));
%     saveas(gcf, fullfile(savepath, [savename '.jpg']));   
%     saveas(gcf, fullfile(savepath, [savename '.epsc']));   
% end    