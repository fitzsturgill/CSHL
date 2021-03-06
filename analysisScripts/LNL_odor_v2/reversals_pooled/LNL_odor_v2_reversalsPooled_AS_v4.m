% this is for noPunish auto reversal blocks...
% normalizing by pre-reversal CSC+ values....
% pertiment to DC_46 and on... I'm carrying through answerLicksROC value...
% LNL_odor_v2_pav_rev  reversal analysis script 9/
%% desktop 
files = {...
    'Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk\DC_46\', 'RE_DC_46.mat';...
    'Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk\DC_47\', 'RE_DC_47.mat';...
%     'Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk\DC_53\', 'RE_DC_53.mat';...
    'Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk\DC_54\', 'RE_DC_54.mat';...
    'Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk\DC_56\', 'RE_DC_56.mat';...
    };


%%
smoothWindow = 3;
%%
for counter = 1:size(files, 1)
    fileName = files{counter, 2}(1:end-4);
    temp = load(fullfile(files{counter, 1}, files{counter, 2}));
    if counter == 1
        RE = temp.RE;
    else
        RE(counter) = temp.RE;
    end
end

%% desktop
savepath = uigetdir;
saveOn = 1;

%% create structure containing all reversals
AR = struct('csPlus', [], 'csMinus', [], 'csPlusReward', [], 'thirdOdor', []); % AR = all reversals
for group = fieldnames(AR)'
    sgroup = group{:};
    for field = fieldnames(RE(1).(sgroup))'
        sfield = field{:};
        for si = 1:length(RE)
            if sum(strcmp(sfield, {'trialsBefore', 'trialsAfter'})) % these are special fields that don't contain before and after data
                continue
            end
            if si == 1
                AR.(sgroup).(sfield).before = RE(si).(sgroup).(sfield).before;
                AR.(sgroup).(sfield).after = RE(si).(sgroup).(sfield).after;
            else
                AR.(sgroup).(sfield).before = expandVertCat(AR.(sgroup).(sfield).before, RE(si).(sgroup).(sfield).before, 'right');
                AR.(sgroup).(sfield).after = expandVertCat(AR.(sgroup).(sfield).after, RE(si).(sgroup).(sfield).after, 'left');
            end
        end
    end
end

firstReversals = cellfun(@(x,y) ~strcmp(x(1:5), y(1:5)), AR.csPlus.filename.after(1:end-1,1), AR.csPlus.filename.after(2:end,1));
firstReversals = [false; firstReversals];
nReversals = size(AR.csPlus.globalTrialNumber.after, 1);
% reversal #s
revNumber = ones(nReversals,1);
mouseNumber = ones(nReversals,1);

thisRev = 1;
thisMouse = 1;
for counter = 2:nReversals
    if strcmp(AR.csPlus.filename.after{counter - 1,1}(1:5), AR.csPlus.filename.after{counter,1}(1:5))
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
comp = {'csLicks', 'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2'}; % comparisons
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

    




% filter reversals according to quality

goodReversals = ...
    ~isnan(trialsToCriterion);% &...
%     auROC.phPeakMean_cs_ch1.before > 0.2 &...
%     auROC.phPeakMean_cs_ch2.before > 0.2;
sortVariable = trialsToCriterion;
[~, sortOrder] = sort(sortVariable);
%
newCsPlus_ch1 = [AR.csMinus.phPeakMean_cs_ch1.before AR.csPlus.phPeakMean_cs_ch1.after];
newCsPlus_ch2 = [AR.csMinus.phPeakMean_cs_ch2.before AR.csPlus.phPeakMean_cs_ch2.after];
newCsPlus_trialNumber = (1:size(newCsPlus_ch1, 2)) - size(AR.csMinus.phPeakMean_cs_ch1.before, 2);
newCsPlus_firstRevTrial = size(AR.csMinus.phPeakMean_cs_ch1.before, 2) + 1;

newCsMinus_ch1 = [AR.csPlus.phPeakMean_cs_ch1.before AR.csMinus.phPeakMean_cs_ch1.after];
newCsMinus_ch2 = [AR.csPlus.phPeakMean_cs_ch2.before AR.csMinus.phPeakMean_cs_ch2.after];
newCsMinus_trialNumber = (1:size(newCsMinus_ch1, 2)) - size(AR.csPlus.phPeakMean_cs_ch1.before, 2);
newCsMinus_firstRevTrial = size(AR.csPlus.phPeakMean_cs_ch1.before, 2) + 1;

alwaysCsPlus_ch1 = [AR.csPlus.phPeakMean_cs_ch1.before AR.csPlus.phPeakMean_cs_ch1.after];
alwaysCsPlus_ch2 = [AR.csPlus.phPeakMean_cs_ch2.before AR.csPlus.phPeakMean_cs_ch2.after];
alwaysCsPlus_trialNumber = (1:size(alwaysCsPlus_ch1, 2)) - size(AR.csPlus.phPeakMean_cs_ch1.before, 2);
alwaysCsPlus_firstRevTrial = size(AR.csPlus.phPeakMean_cs_ch1.before, 2) + 1;

odor3_ch1 = [AR.thirdOdor.phPeakMean_cs_ch1.before AR.thirdOdor.phPeakMean_cs_ch1.after];
odor3_ch2 = [AR.thirdOdor.phPeakMean_cs_ch2.before AR.thirdOdor.phPeakMean_cs_ch2.after];
odor3_licks = [AR.thirdOdor.csLicks.before AR.thirdOdor.csLicks.after];
odor3_trialNumber = (1:size(odor3_ch1, 2)) - size(AR.thirdOdor.phPeakMean_cs_ch1.before, 2);

oldCsPlus_trialNumber = -(size(AR.csPlus.phPeakMean_cs_ch1.before, 2) - 1) : 0;
oldCsMinus_trialNumber = -(size(AR.csMinus.phPeakMean_cs_ch1.before, 2) - 1) : 0;

% compile data

newCsPlus_roc = [AR.csMinus.csLicksROC.before AR.csPlus.csLicksROC.after];
newCsMinus_roc = [AR.csPlus.csLicksROC.before AR.csMinus.csLicksROC.after];
alwaysCsPlus_roc = [AR.csPlus.csLicksROC.before AR.csPlus.csLicksROC.after];

newCsPlus_licks = [AR.csMinus.csLicks.before AR.csPlus.csLicks.after];
newCsMinus_licks = [AR.csPlus.csLicks.before AR.csMinus.csLicks.after];
alwaysCsPlus_licks = [AR.csPlus.csLicks.before AR.csPlus.csLicks.after];

newCsPlus_ch1_norm = smoothdata(newCsPlus_ch1, 2, 'movmean', smoothWindow, 'omitnan');
newCsPlus_ch2_norm = smoothdata(newCsPlus_ch2, 2, 'movmean', smoothWindow, 'omitnan');
newCsPlus_licks_norm = smoothdata(newCsPlus_licks, 2, 'movmean', smoothWindow, 'omitnan');

newCsMinus_ch1_norm = smoothdata(newCsMinus_ch1, 2, 'movmean', smoothWindow, 'omitnan');
newCsMinus_ch2_norm = smoothdata(newCsMinus_ch2, 2, 'movmean', smoothWindow, 'omitnan');
newCsMinus_licks_norm = smoothdata(newCsMinus_licks, 2, 'movmean', smoothWindow, 'omitnan');

alwaysCsPlus_ch1_norm = smoothdata(alwaysCsPlus_ch1, 2, 'movmean', smoothWindow, 'omitnan');
alwaysCsPlus_ch2_norm = smoothdata(alwaysCsPlus_ch2, 2, 'movmean', smoothWindow, 'omitnan');
alwaysCsPlus_licks_norm = smoothdata(alwaysCsPlus_licks, 2, 'movmean', smoothWindow, 'omitnan');

% alwaysCsPlus_ch3 = 
% normalize by f% of pre reversal values, smooth, find first and last common points across reversals

f = 0.8;
trialsBack = 20;
%
% IT SEEMS LIKE WHAT WORKS BEST IS JUST TO NORMALIZE USING THE 80% OR SO
% OF PRE-REVERSAL CS+ VALUES

% normalize by csPlus before reversal
% normVector_ch1 = percentile(newCsMinus_ch1_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), f, 2);
% normVector_ch2 = percentile(newCsMinus_ch2_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), f, 2);
% normVector_licks = percentile(newCsMinus_licks_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), f, 2);
normVector_ch1 = nanmean(newCsMinus_ch1_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), 2);
normVector_ch2 = nanmean(newCsMinus_ch2_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), 2);
normVector_licks = nanmedian(newCsMinus_licks_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), 2);

% alternatively, try normalizing by Odor3 cue response
% normVector_ch1 = nanmedian(odor3_ch1, 2);
% normVector_ch2 = nanmean(odor3_ch2, 2);
% normVector_licks = nanmean(odor3_licks, 2);


% only normalization 
newCsPlus_ch1_norm = newCsPlus_ch1_norm ./ normVector_ch1; % singleton expansion by default
newCsPlus_ch2_norm = newCsPlus_ch2_norm ./ normVector_ch2;
newCsPlus_licks_norm = newCsPlus_licks_norm ./ normVector_licks;

newCsMinus_ch1_norm = newCsMinus_ch1_norm ./ normVector_ch1;
newCsMinus_ch2_norm = newCsMinus_ch2_norm ./ normVector_ch2;
newCsMinus_licks_norm = newCsMinus_licks_norm ./ normVector_licks;

alwaysCsPlus_ch1_norm = alwaysCsPlus_ch1_norm ./ normVector_ch1;
alwaysCsPlus_ch2_norm = alwaysCsPlus_ch2_norm ./ normVector_ch2;
alwaysCsPlus_licks_norm = alwaysCsPlus_licks_norm ./ normVector_licks;

odor3_ch1_norm = odor3_ch1 ./ normVector_ch1;
odor3_ch2_norm = odor3_ch2 ./ normVector_ch2;
odor3_licks_norm = odor3_licks ./ normVector_licks;

%% make wrap-around data arrays
%%


% % combined tiled plot
% sb = [3 3];
% savename1 = 'newCsPlus_combined_tiled';
% 
% h1 = ensureFigure(savename1, 1);
% mcLandscapeFigSetup(h1);
% for counter = 1:9    
%     a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); 
%     counter = counter + 18;
%     plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch1(counter, :), 'movmean', smoothWindow, 'omitnan'), '-g'); hold on
%     plot(oldCsPlus_trialNumber, smoothdata(AR.csPlus.phPeakMean_cs_ch1.before(counter, :), 'movmean', smoothWindow, 'omitnan'), '--g');
%     plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch2(counter, :), 'movmean', smoothWindow, 'omitnan'), '-r'); hold on
%     plot(oldCsPlus_trialNumber, smoothdata(AR.csPlus.phPeakMean_cs_ch2.before(counter, :), 'movmean', smoothWindow, 'omitnan'), '--r');
%     textBox(num2str(counter));
% end

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
%% annotate good reversals and plot individual reversals
tilePos = [961.0000  955.4000  876.0000  648.0000];
reversalsColor = repmat([0 1 0], nReversals, 1);
reversalsColor = bsxfun(@times, reversalsColor, goodReversals);
sb = repmat(ceil(sqrt(nReversals)), 1, 2);
savename1 = 'newCsPlus_ch1_tiled';
h1 = ensureFigure(savename1, 1);
savename2 = 'newCsPlus_ch2_tiled';
h2 = ensureFigure(savename2, 1);
savename3 = 'newCsPlus_licks_tiled';
h3 = ensureFigure(savename3, 1);
for counter = 1:nReversals
    a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); 
    plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch1(counter, :), 'movmean', 3, 'omitnan'), 'Parent', a1, 'Color', reversalsColor(counter,:));
    a2 = subplot(sb(1),sb(2),counter, 'Parent', h2); 
    plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch2(counter, :), 'movmean', 3, 'omitnan'), 'Parent', a2, 'Color', reversalsColor(counter,:));    
    a3 = subplot(sb(1),sb(2),counter, 'Parent', h3); 
    plot(newCsPlus_trialNumber, smoothdata(newCsPlus_licks(counter, :), 'movmean', 3, 'omitnan'), 'Parent', a3, 'Color', reversalsColor(counter,:));        
end
% set([h1 h2 h3], 'Position', tilePos);

sb = repmat(ceil(sqrt(nReversals)), 1, 2);
savename1 = 'newCsMinus_ch1_tiled';
h1 = ensureFigure(savename1, 1);
savename2 = 'newCsMinus_ch2_tiled';
h2 = ensureFigure(savename2, 1);
for counter = 1:nReversals
    a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); 
    plot(newCsMinus_trialNumber, smoothdata(newCsMinus_ch1(counter, :), 'movmean', 3, 'omitnan'), 'Parent', a1, 'Color', reversalsColor(counter,:));
    a2 = subplot(sb(1),sb(2),counter, 'Parent', h2); 
    plot(newCsMinus_trialNumber, smoothdata(newCsMinus_ch2(counter, :), 'movmean', 3, 'omitnan'), 'Parent', a2, 'Color', reversalsColor(counter,:));    
end
% set([h1 h2], 'Position', tilePos);

sb = repmat(ceil(sqrt(nReversals)), 1, 2);
savename1 = 'alwaysCsPlus_ch1_tiled';
h1 = ensureFigure(savename1, 1);
savename2 = 'alwaysCsPlus_ch2_tiled';
h2 = ensureFigure(savename2, 1);
for counter = 1:nReversals
    a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); 
    plot(alwaysCsPlus_trialNumber, smoothdata(alwaysCsPlus_ch1(counter, :), 'movmean', 1, 'omitnan'), 'Parent', a1, 'Color', reversalsColor(counter,:));
    a2 = subplot(sb(1),sb(2),counter, 'Parent', h2); 
    plot(alwaysCsPlus_trialNumber, smoothdata(alwaysCsPlus_ch2(counter, :), 'movmean', 1, 'omitnan'), 'Parent', a2, 'Color', reversalsColor(counter,:));    
end
% set([h1 h2], 'Position', tilePos);

sb = repmat(ceil(sqrt(nReversals)), 1, 2);
savename1 = 'before_ch1_tiled';
h1 = ensureFigure(savename1, 1);
savename2 = 'before_ch2_tiled';
h2 = ensureFigure(savename2, 1);
for counter = 1:nReversals
    a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); hold(a1, 'on');
    plot(oldCsPlus_trialNumber, smoothdata(AR.csPlus.phPeakMean_cs_ch1.before(counter, :), 'movmean', smoothWindow, 'omitnan'), 'Parent', a1, 'Color', reversalsColor(counter,:)); 
    plot(oldCsMinus_trialNumber, smoothdata(AR.csMinus.phPeakMean_cs_ch1.before(counter, :), 'movmean', smoothWindow, 'omitnan'), 'Parent', a1, 'LineStyle', '--', 'Color', reversalsColor(counter,:)); 
    a2 = subplot(sb(1),sb(2),counter, 'Parent', h2); hold(a2, 'on');
    plot(oldCsPlus_trialNumber, smoothdata(AR.csPlus.phPeakMean_cs_ch2.before(counter, :), 'movmean', smoothWindow, 'omitnan'), 'Parent', a2, 'Color', reversalsColor(counter,:)); 
    plot(oldCsMinus_trialNumber, smoothdata(AR.csMinus.phPeakMean_cs_ch2.before(counter, :), 'movmean', smoothWindow, 'omitnan'), 'Parent', a2, 'LineStyle', '--', 'Color', reversalsColor(counter,:)); 
end
% set([h1 h2], 'Position', tilePos);

%% images

% newCsPlus
saveName = 'newCsPlus_image';
ensureFigure(saveName, 1);
subplot(2,2,1);
clim = [-2 2];
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

%% normalized
common = sum(~isnan(newCsPlus_ch1_norm)) > 3;
common_odor3 = (-9 <= odor3_trialNumber) & (odor3_trialNumber <= 9);

xlim = [-30 30];
ylim = [-1 2];

savename = 'reversals_newCsPlus';
ensureFigure(savename, 1);
hla = zeros(1,6);
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_ch2_norm(goodReversals, common)), nanSEM(newCsPlus_ch2_norm(goodReversals, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_ch1_norm(goodReversals, common)), nanSEM(newCsPlus_ch1_norm(goodReversals, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmedian(newCsPlus_licks_norm(goodReversals, common)), zeros(1, sum(common)),...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch1_norm(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch2_norm(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmedian(odor3_licks_norm(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
    set(gca, 'XLim', xlim, 'YLim', ylim);
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
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmedian(newCsMinus_licks_norm(goodReversals, newCsMinus_common)), zeros(1, sum(common)),...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch1_norm(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch2_norm(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmedian(odor3_licks_norm(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
    set(gca, 'XLim', xlim, 'YLim', ylim);    
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
o3xdata = (0:(size(AR.thirdOdor.csLicks.before, 2) + size(AR.thirdOdor.csLicks.after, 2) - 1)) - size(AR.thirdOdor.csLicks.before, 2) - 1; 
o3licks = [AR.thirdOdor.csLicks.before AR.thirdOdor.csLicks.after];
o3_ch1 = [AR.thirdOdor.phPeakMean_cs_ch1.before AR.thirdOdor.phPeakMean_cs_ch1.after];
o3_ch2 = [AR.thirdOdor.phPeakMean_cs_ch2.before AR.thirdOdor.phPeakMean_cs_ch2.after];
hla = zeros(1,6);
[hl, hp] = boundedline(o3xdata, nanmean(odor3_licks_norm(goodReversals, :)), nanSEM(odor3_licks_norm(goodReversals, :))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(o3xdata, nanmean(odor3_ch1_norm(goodReversals, :)), nanSEM(odor3_ch1_norm(goodReversals, :))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(o3xdata, nanmean(odor3_ch2_norm(goodReversals, :)), nanSEM(odor3_ch2_norm(goodReversals, :))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch1_norm(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch2_norm(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmedian(odor3_licks_norm(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);    
    set(gca, 'XLim', [-10 10]);%, 'YLim', [-1 2]);
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks'...
                '\bf\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\bf\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\bf\color[rgb]{0.5,0.5,0.5}Odor 3'},...
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
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmedian(alwaysCsPlus_licks_norm(goodReversals, common)), zeros(1, sum(common)),...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch1_norm(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3_ch2_norm(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmedian(odor3_licks_norm(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
    set(gca, 'XLim', xlim, 'YLim', ylim);    
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


%% find good example reversals for Ach, Dop, and licking
sb = repmat(ceil(sqrt(nReversals)), 1, 2);
savename1 = 'newCsPlus_all_tiled';
h1 = ensureFigure(savename1, 1);
savename2 = 'newCsMinus_all_tiled';
h2 = ensureFigure(savename2, 1);
for counter = 1:nReversals
    a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); hold(a1, 'on');
    plot(newCsPlus_trialNumber, newCsPlus_ch1_norm(counter, :), 'g', 'Parent', a1); 
    plot(newCsPlus_trialNumber, newCsPlus_ch2_norm(counter, :), 'r', 'Parent', a1);
    plot(newCsPlus_trialNumber, newCsPlus_licks_norm(counter, :), 'k', 'Parent', a1);    
    a2 = subplot(sb(1),sb(2),counter, 'Parent', h2); hold(a2, 'on');
    plot(newCsMinus_trialNumber, newCsMinus_ch1_norm(counter, :), 'g', 'Parent', a2); 
    plot(newCsMinus_trialNumber, newCsMinus_ch2_norm(counter, :), 'r', 'Parent', a2);
    plot(newCsMinus_trialNumber, newCsMinus_licks_norm(counter, :), 'k', 'Parent', a2);        
end
if saveOn
    saveas(h1, fullfile(savepath, [savename1 '.fig']));
    saveas(h1, fullfile(savepath, [savename1 '.jpg']));   
    saveas(h1, fullfile(savepath, [savename1 '.epsc']));   
    saveas(h2, fullfile(savepath, [savename2 '.fig']));
    saveas(h2, fullfile(savepath, [savename2 '.jpg']));   
    saveas(h2, fullfile(savepath, [savename2 '.epsc']));   
end    

%% plot good example reversals for both conditions  (CS- -> CS+,  CS+ -> CS-) 
ex1 = 11;
ex2 = 2; %3;
acolor = [0.6680,0.2148,0.8359];
dcolor = [0.9258,0.4883,0.1914];
lcolor = [0.3 0.3 0.3];
lwidth = 2;
savename = 'newCsPlus_example';
smoothFactor = 5;
smoothMethod = 'movmean';
ensureFigure(savename, 1); axes; hold on;
plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch1_norm(ex1, :), smoothMethod, smoothFactor), 'Color', acolor, 'LineWidth', lwidth);
plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch2_norm(ex1, :), smoothMethod, smoothFactor), 'Color', dcolor, 'LineWidth', lwidth);
plot(newCsPlus_trialNumber, smoothdata(newCsPlus_licks_norm(ex1, :), smoothMethod, smoothFactor), 'Color', lcolor, 'LineWidth', lwidth);
set(gca, 'XLim', [-40 40], 'XTick', [], 'YLim', [-1 1.5]);
addOrginLines;
formatFigurePoster([5.5 2], '', 20);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

savename = 'newCsMinus_example';
ensureFigure(savename, 1); axes; hold on;
plot(newCsMinus_trialNumber, smoothdata(newCsMinus_ch1_norm(ex2, :), smoothMethod, smoothFactor), 'Color', acolor, 'LineWidth', lwidth);
plot(newCsMinus_trialNumber, smoothdata(newCsMinus_ch2_norm(ex2, :), smoothMethod, smoothFactor), 'Color', dcolor, 'LineWidth', lwidth);
plot(newCsMinus_trialNumber, smoothdata(newCsMinus_licks_norm(ex2, :), smoothMethod, smoothFactor), 'Color', lcolor, 'LineWidth', lwidth);
set(gca, 'XLim', [-40 40], 'XTick', [], 'YLim', [-1 1.5]);
addOrginLines;
formatFigurePoster([5.5 2], '', 20);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    