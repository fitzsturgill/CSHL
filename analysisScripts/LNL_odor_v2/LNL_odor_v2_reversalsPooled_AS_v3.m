% normalizing by pre-reversal CSC+ values....
% LNL_odor_v2_pav_rev  reversal analysis script 9/
%% desktop 
files = {...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_17\', 'RE_DC_17.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_20\', 'RE_DC_20.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_36\', 'RE_DC_36.mat';...    
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_37\', 'RE_DC_37.mat';...    
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_40\', 'RE_DC_40.mat';... % first reversal currently excluded on DC_40
    };

%% laptop

% files = {...
%     'C:\Fitz_Data\SummaryAnalyses\LNL_Analysis_new\DC_17\', 'RE_DC_17.mat';...
%     'C:\Fitz_Data\SummaryAnalyses\LNL_Analysis_new\DC_20\', 'RE_DC_20.mat';...
%     };

% files = {...
%     'C:\Fitz_Data\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_17\', 'RE_DC_17.mat';...
%     'C:\Fitz_Data\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_20\', 'RE_DC_20.mat';...
%     'C:\Fitz_Data\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_36\', 'RE_DC_36.mat';...    
% %     'C:\Fitz_Data\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_37\', 'RE_DC_37.mat';...    
%     'C:\Fitz_Data\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_40\', 'RE_DC_40.mat';... % first reversal currently excluded on DC_40
%     };

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
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SFN_2017\Reversals';
savepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\Reversals_Pooled';
saveOn = 1;






%% create structure containing all reversals
AR = struct('csPlus', [], 'csMinus', [], 'csPlusReward', []); % AR = all reversals
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
[~, sortOrder] = sort(revNumber);
    




%% quality control- calculate auROC and dPrime for relevent comparisons
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
    auROC.csLicks.before > 0 &...
    auROC.csLicks.after > -0.7 &...
    auROC.phPeakMean_cs_ch1.before > 0.5 &...
    auROC.phPeakMean_cs_ch2.before > 0.5;

% goodReversals = auROC.csLicks.acq > 0.5 & auROC.phPeakMean_cs_ch1.before > 0.4 & auROC.phPeakMean_cs_ch2.before > 0.4;


%%
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

oldCsPlus_trialNumber = -(size(AR.csPlus.phPeakMean_cs_ch1.before, 2) - 1) : 0;
oldCsMinus_trialNumber = -(size(AR.csMinus.phPeakMean_cs_ch1.before, 2) - 1) : 0;

%% compile data
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

% normalize by f% of pre reversal values, smooth, find first and last common points across reversals

f = 0.8;
trialsBack = 50;

% normalize by csPlus before reversal
normVector_ch1 = percentile(smoothdata(newCsMinus_ch1_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), 2, 'movmean', 3, 'omitnan'), f, 2);
normVector_ch2 = percentile(smoothdata(newCsMinus_ch2_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), 2, 'movmean', 3, 'omitnan'), f, 2);
normVector_licks = percentile(smoothdata(newCsMinus_licks_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), 2, 'movmean', 3, 'omitnan'), f, 2);

% normalize by csPlus after reversal
% normVector_ch1 = percentile(newCsPlus_ch1_norm(:,newCsPlus_firstRevTrial:end), f, 2);
% normVector_ch2 = percentile(newCsPlus_ch2_norm(:,newCsPlus_firstRevTrial:end), f, 2);
% normVector_licks = percentile(newCsPlus_licks_norm(:,newCsPlus_firstRevTrial:end), f, 2);

% no normalization
% normVector_ch1 = ones(nReversals, 1);
% normVector_ch2 = ones(nReversals, 1);
% normVector_licks = ones(nReversals, 1);

newCsPlus_ch1_norm = bsxfun(@rdivide, newCsPlus_ch1_norm, normVector_ch1);
newCsPlus_ch2_norm = bsxfun(@rdivide, newCsPlus_ch2_norm, normVector_ch2);
newCsPlus_licks_norm = bsxfun(@rdivide, newCsPlus_licks_norm, normVector_licks);

newCsMinus_ch1_norm = bsxfun(@rdivide, newCsMinus_ch1_norm, normVector_ch1);
newCsMinus_ch2_norm = bsxfun(@rdivide, newCsMinus_ch2_norm, normVector_ch2);
newCsMinus_licks_norm = bsxfun(@rdivide, newCsMinus_licks_norm, normVector_licks);

alwaysCsPlus_ch1_norm = bsxfun(@rdivide, alwaysCsPlus_ch1_norm, normVector_ch1);
alwaysCsPlus_ch2_norm = bsxfun(@rdivide, alwaysCsPlus_ch2_norm, normVector_ch2);
alwaysCsPlus_licks_norm = bsxfun(@rdivide, alwaysCsPlus_licks_norm, normVector_licks);

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
set([h1 h2 h3], 'Position', tilePos);

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
set([h1 h2], 'Position', tilePos);

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
set([h1 h2], 'Position', tilePos);

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
set([h1 h2], 'Position', tilePos);






%% images

% newCsPlus
ensureFigure('newCsPlus_image', 1);
subplot(2,2,1);
clim = [-10 10];
xlim = [min(newCsPlus_trialNumber), max(newCsPlus_trialNumber)];
imagesc('XData', xlim, 'CData', newCsPlus_ch1_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('licks');  set(gca, 'CLim', [-5 5])
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
subplot(2,2,2);
imagesc('XData', xlim, 'CData', newCsPlus_ch2_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Ach');  set(gca, 'CLim', [-5 5])
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
subplot(2,2,3);
imagesc('XData', xlim, 'CData', newCsPlus_licks_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on;title('Dop');  set(gca, 'CLim', [-5 5])
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');

% newCsMinus
ensureFigure('newCsMinus_image', 1);
subplot(2,2,1);
xlim = [min(newCsMinus_trialNumber), max(newCsMinus_trialNumber)];
imagesc('XData', xlim, 'CData', newCsMinus_ch1_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('licks');  set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
subplot(2,2,2);
imagesc('XData', xlim, 'CData', newCsMinus_ch2_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Ach');  set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
subplot(2,2,3);
imagesc('XData', xlim, 'CData', newCsMinus_licks_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Dop');  set(gca, 'CLim', clim)
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');



% alwaysCsPlus 
ensureFigure('alwaysCsPlus_image', 1);
subplot(2,2,1);
xlim = [min(alwaysCsPlus_trialNumber), max(alwaysCsPlus_trialNumber)];
imagesc('XData', xlim, 'CData', alwaysCsPlus_ch1_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('licks');  set(gca, 'CLim', [-1.5 1.5])
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
subplot(2,2,2);
imagesc('XData', xlim, 'CData', alwaysCsPlus_ch2_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Ach');  set(gca, 'CLim', [-1 1])
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');
subplot(2,2,3);
imagesc('XData', xlim, 'CData', alwaysCsPlus_licks_norm(sortOrder, :)); set(gca, 'XLim', xlim); hold on; title('Dop');  set(gca, 'CLim', [-1.5 1.5])
scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:nReversals, [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled');

%%
common = sum(~isnan(newCsPlus_ch1_norm)) > 3;


savename = 'reversals_newCsPlus';
ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_ch2_norm(goodReversals, common)), nanSEM(newCsPlus_ch2_norm(goodReversals, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_ch1_norm(goodReversals, common)), nanSEM(newCsPlus_ch1_norm(goodReversals, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_licks_norm(goodReversals, common)), nanSEM(newCsPlus_licks_norm(goodReversals, common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
set(hla, 'LineWidth', 2);
set(gca, 'XLim', [-50 80]);%, 'YLim', [-1 1]);
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', ...
        '\bf\color[rgb]{0.5,0.5,0.5}Licks'}, 'Location', 'southeast', 'FontSize', 18, 'Interpreter', 'tex', 'Box', 'off');

    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');
    
    formatFigurePoster([5.5 4], '', 20);
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
hla = zeros(1,3);
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmean(newCsMinus_ch2_norm(goodReversals, newCsMinus_common)), nanSEM(newCsMinus_ch2_norm(goodReversals, newCsMinus_common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmean(newCsMinus_ch1_norm(goodReversals, newCsMinus_common)), nanSEM(newCsMinus_ch1_norm(goodReversals, newCsMinus_common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmean(newCsMinus_licks_norm(goodReversals, newCsMinus_common)), nanSEM(newCsMinus_licks_norm(goodReversals, newCsMinus_common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
set(hla, 'LineWidth', 2);
% set(gca, 'XLim', [-50 80]);%, 'YLim', [-1 1]);
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', ...
        '\bf\color[rgb]{0.5,0.5,0.5}Licks'}, 'Location', 'southeast', 'FontSize', 18, 'Interpreter', 'tex', 'Box', 'off');

    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');
    
    formatFigurePoster([5.5 4], '', 20);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    


%%

% common = find(mean(alwaysCsPlus_ch1_norm));% applying mean will reveal NaNs
common = sum(~isnan(alwaysCsPlus_ch1_norm)) > 3;% applying mean will reveal NaNs



savename = 'reversals_alwaysCsPlus';
ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus_ch2_norm(goodReversals, common)), nanSEM(alwaysCsPlus_ch2_norm(goodReversals, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus_ch1_norm(goodReversals, common)), nanSEM(alwaysCsPlus_ch1_norm(goodReversals, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus_licks_norm(goodReversals, common)), nanSEM(alwaysCsPlus_licks_norm(goodReversals, common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
set(hla, 'LineWidth', 2);
% set(gca, 'XLim', [-40 40]);%, 'YLim', [-1 1.6]);
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', ...
        '\bf\color[rgb]{0.5,0.5,0.5}Licks'}, 'Location', 'southeast', 'FontSize', 18, 'Interpreter', 'tex', 'Box', 'off');

    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');
    
    formatFigurePoster([5.5 4], '', 20);
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