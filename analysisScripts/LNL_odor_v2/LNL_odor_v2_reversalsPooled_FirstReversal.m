%% desktop 
files = {...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_17\', 'RE_DC_17.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_18\', 'RE_DC_18.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_20\', 'RE_DC_20.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_25\', 'RE_DC_25.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_27\', 'RE_DC_27.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_36\', 'RE_DC_36.mat';...
    'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal\DC_40\', 'RE_DC_40.mat';...
    };


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

%% 
%% desktop
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SFN_2017\Reversals';
savepath = 'Z:\SummaryAnalyses\LNL_odor_v2_BaselineTrialByTrial_firstReversal';
saveOn = 1;

%% create structure containing all reversals
AR = struct('csPlus', [], 'csMinus', [], 'csPlusReward', [], 'allTrials', []); % AR = all reversals
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
nReversals = size(AR.csPlus.globalTrialNumber.after, 1);
if saveOn
    save(fullfile(savepath, 'AR.mat'), 'AR');
    disp(['*** Saved: ' fullfile(savepath, 'AR.mat')]);
end

%% figure out indices of first reversals, and valid experiments for photometry data and compile data
firstReversals = cellfun(@(x,y) ~strcmp(x(1:5), y(1:5)), AR.csPlus.filename.after(1:end-1,1), AR.csPlus.filename.after(2:end,1));
firstReversals = [false; firstReversals];

% Ch1List = {'DC_17', 'DC_20', 'DC_36', 'DC_37', 'DC_40'};
% Ch2List = {'DC_17', 'DC_20', 'DC_25', 'DC_27', 'DC_36', 'DC_37', 'DC_40'}; % includes VTA fiber of dual 6f mice

ch1Reversals = cellfun(@(x) any(strcmp(x(1:5), {'DC_17', 'DC_20', 'DC_36', 'DC_37', 'DC_40'})), AR.csPlus.filename.after(:,1));
ch2Reversals = cellfun(@(x) any(strcmp(x(1:5), {'DC_17', 'DC_20', 'DC_25', 'DC_27', 'DC_36', 'DC_37', 'DC_40'})), AR.csPlus.filename.after(:,1));

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

newCsPlus_cueLicks = [AR.csMinus.csLicks.before AR.csPlus.csLicks.after];
newCsMinus_cueLicks = [AR.csPlus.csLicks.before AR.csMinus.csLicks.after];
alwaysCsPlus_cueLicks = [AR.csPlus.csLicks.before AR.csPlus.csLicks.after];

%% newCsPlus and newCsMinus - normalize by 90% of pre reversal CS+ values, smooth, find first and last common points across reversals
f = 0.9;
newCsPlus_ch1_norm = smoothdata(newCsPlus_ch1, 2, 'movmean', 3, 'omitnan');
newCsPlus_ch2_norm = smoothdata(newCsPlus_ch2, 2, 'movmean', 3, 'omitnan');
newCsPlus_cueLicks_norm = smoothdata(newCsPlus_cueLicks, 2, 'movmean', 3, 'omitnan');



newCsPlus_ch1_norm = bsxfun(@rdivide, newCsPlus_ch1_norm, percentile(newCsMinus_ch1_norm(:,1:newCsMinus_firstRevTrial - 1), f, 2));
newCsPlus_ch2_norm = bsxfun(@rdivide, newCsPlus_ch2_norm, percentile(newCsMinus_ch2_norm(:,1:newCsMinus_firstRevTrial - 1), f, 2));


newCsMinus_ch1_norm = smoothdata(newCsMinus_ch1, 2, 'movmean', 3, 'omitnan');
newCsMinus_ch2_norm = smoothdata(newCsMinus_ch2, 2, 'movmean', 3, 'omitnan');
newCsMinus_cueLicks_norm = smoothdata(newCsMinus_cueLicks, 2, 'movmean', 3, 'omitnan');

newCsMinus_ch1_norm = bsxfun(@rdivide, newCsMinus_ch1_norm, percentile(newCsMinus_ch1_norm(:,1:newCsMinus_firstRevTrial - 1), f, 2));
newCsMinus_ch2_norm = bsxfun(@rdivide, newCsMinus_ch2_norm, percentile(newCsMinus_ch2_norm(:,1:newCsMinus_firstRevTrial - 1), f, 2));
    

% find first pseudo-common trials across reversals, plot antic lick rate for first reversal vs. subsequent reversals
% common_newCsPlus = sum(~isnan(newCsPlus_ch1_norm)) > 3;
% common_newCsMinus = sum(~isnan(newCsMinus_ch1_norm)) > 3;

common_newCsPlus = false([1 size(newCsPlus_ch1_norm, 2)]);
common_newCsPlus(newCsPlus_firstRevTrial - 50 : newCsPlus_firstRevTrial + 50) = true;
common_newCsMinus = false([1 size(newCsMinus_ch1_norm, 2)]);
common_newCsMinus(newCsPlus_firstRevTrial - 50 : newCsMinus_firstRevTrial + 50) = true;
%%
trialwindow = [-50 50];
savename = 'first_reversal_vs_others';
ensureFigure(savename, 1); hold on;
subplot(3,2,1);
hla = zeros(1,2);
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsPlus), nanmean(newCsPlus_cueLicks_norm(firstReversals, common_newCsPlus)),...
    nanSEM(newCsPlus_cueLicks_norm(firstReversals, common_newCsPlus))', 'cmap', [0.5 0.5 0.5]);
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsPlus), nanmean(newCsPlus_cueLicks_norm(~firstReversals, common_newCsPlus)), nanSEM(newCsPlus_cueLicks_norm(~firstReversals, common_newCsPlus))',...
    'cmap', [0 0 0]);
hla(2) = hl;
set(hla, 'LineWidth', 2);
set(hla(1), 'LineStyle', '--'); 
set(gca, 'XLim', trialwindow);
h  = addOrginLines;
ylabel('Cue licks/s');
title('New Cs+');

% licks, new Cs-
subplot(3,2,2);
hla = zeros(1,2);
[hl, hp] = boundedline(newCsMinus_trialNumber(common_newCsMinus), nanmean(newCsMinus_cueLicks_norm(firstReversals, common_newCsMinus)),...
    nanSEM(newCsMinus_cueLicks_norm(firstReversals, common_newCsMinus))', 'cmap', [0.5 0.5 0.5]);
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(common_newCsMinus), nanmean(newCsMinus_cueLicks_norm(~firstReversals, common_newCsMinus)), nanSEM(newCsMinus_cueLicks_norm(~firstReversals, common_newCsMinus))',...
    'cmap', [0 0 0]);
hla(2) = hl;
set(hla, 'LineWidth', 2);
set(hla(1), 'LineStyle', '--'); 
set(gca, 'XLim', trialwindow);
h  = addOrginLines;
title('New Cs-');

subplot(3,2,3);
hla = zeros(1,2);
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsPlus), nanmean(newCsPlus_ch1_norm(firstReversals & ch1Reversals, common_newCsPlus)), nanSEM(newCsPlus_ch1_norm(firstReversals & ch1Reversals, common_newCsPlus))',...
    'cmap', [171 55 214]/256 * 0.5);
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsPlus), nanmean(newCsPlus_ch1_norm(~firstReversals & ch1Reversals, common_newCsPlus)), nanSEM(newCsPlus_ch1_norm(~firstReversals & ch1Reversals, common_newCsPlus))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
set(hla, 'LineWidth', 2);
set(hla(1), 'LineStyle', '--'); 
set(gca, 'XLim', trialwindow);
h  = addOrginLines;
ylabel('BF ACh.');

subplot(3,2,4);
hla = zeros(1,2);
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsMinus), nanmean(newCsMinus_ch1_norm(firstReversals & ch1Reversals, common_newCsMinus)), nanSEM(newCsMinus_ch1_norm(firstReversals & ch1Reversals, common_newCsMinus))',...
    'cmap', [171 55 214]/256 * 0.5);
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsMinus), nanmean(newCsMinus_ch1_norm(~firstReversals & ch1Reversals, common_newCsMinus)), nanSEM(newCsMinus_ch1_norm(~firstReversals & ch1Reversals, common_newCsMinus))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
set(hla, 'LineWidth', 2);
set(hla(1), 'LineStyle', '--'); 
set(gca, 'XLim', trialwindow);
h  = addOrginLines;

subplot(3,2,5);
hla = zeros(1,2);
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsPlus), nanmean(newCsPlus_ch2_norm(firstReversals & ch2Reversals, common_newCsPlus)), nanSEM(newCsPlus_ch2_norm(firstReversals & ch2Reversals, common_newCsPlus))',...
    'cmap', [237 125 49]/256 * 0.5);
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsPlus), nanmean(newCsPlus_ch2_norm(~firstReversals & ch2Reversals, common_newCsPlus)), nanSEM(newCsPlus_ch2_norm(~firstReversals & ch2Reversals, common_newCsPlus))',...
    'cmap', [237 125 49]/256);
hla(2) = hl;
set(hla, 'LineWidth', 2);
set(hla(1), 'LineStyle', '--'); 
set(gca, 'XLim', trialwindow);
h  = addOrginLines;
ylabel('VTA Dop.'); xlabel('Odor Presentations from Reversal');

subplot(3,2,6);
hla = zeros(1,2);
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsPlus), nanmean(newCsMinus_ch2_norm(firstReversals & ch2Reversals, common_newCsPlus)), nanSEM(newCsMinus_ch2_norm(firstReversals & ch2Reversals, common_newCsPlus))',...
    'cmap', [237 125 49]/256 * 0.5);
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common_newCsPlus), nanmean(newCsMinus_ch2_norm(~firstReversals & ch2Reversals, common_newCsPlus)), nanSEM(newCsMinus_ch2_norm(~firstReversals & ch2Reversals, common_newCsPlus))',...
    'cmap', [237 125 49]/256);
hla(2) = hl;
set(hla, 'LineWidth', 2);
set(hla(1), 'LineStyle', '--'); 
set(gca, 'XLim', trialwindow);
h  = addOrginLines; xlabel('Odor Presentations from Reversal');
    

formatFigurePoster([10 8], '', 14);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    