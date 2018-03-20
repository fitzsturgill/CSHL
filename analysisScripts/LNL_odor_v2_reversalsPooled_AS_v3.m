% normalizing by pre-reversal CSC+ values....
% LNL_odor_v2_pav_rev  reversal analysis script 9/
%% desktop 
files = {...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_17\', 'RE_DC_17.mat';...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_20\', 'RE_DC_20.mat';...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_35\', 'RE_DC_35.mat';...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_36\', 'RE_DC_36.mat';...    
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_37\', 'RE_DC_37.mat';...    
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_40\', 'RE_DC_40.mat';... % first reversal currently excluded on DC_40
    };

%% laptop

% files = {...
%     'C:\Fitz_Data\SummaryAnalyses\LNL_Analysis_new\DC_17\', 'RE_DC_17.mat';...
%     'C:\Fitz_Data\SummaryAnalyses\LNL_Analysis_new\DC_20\', 'RE_DC_20.mat';...
%     };

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
nReversals = size(AR.csPlus.globalTrialNumber.after, 1);
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

sb = repmat(ceil(sqrt(nReversals)), 1, 2);
savename1 = 'newCsPlus_ch1_tiled';
h1 = ensureFigure(savename1, 1);
savename2 = 'newCsPlus_ch2_tiled';
h2 = ensureFigure(savename2, 1);
for counter = 1:nReversals
    a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); 
    plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch1(counter, :), 'movmean', 3, 'omitnan'), 'Parent', a1);
    a2 = subplot(sb(1),sb(2),counter, 'Parent', h2); 
    plot(newCsPlus_trialNumber, smoothdata(newCsPlus_ch2(counter, :), 'movmean', 3, 'omitnan'), 'Parent', a2);    
end
sb = repmat(ceil(sqrt(nReversals)), 1, 2);
savename1 = 'newCsMinus_ch1_tiled';
h1 = ensureFigure(savename1, 1);
savename2 = 'newCsMinus_ch2_tiled';
h2 = ensureFigure(savename2, 1);
for counter = 1:nReversals
    a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); 
    plot(newCsMinus_trialNumber, smoothdata(newCsMinus_ch1(counter, :), 'movmean', 3, 'omitnan'), 'Parent', a1);
    a2 = subplot(sb(1),sb(2),counter, 'Parent', h2); 
    plot(newCsMinus_trialNumber, smoothdata(newCsMinus_ch2(counter, :), 'movmean', 3, 'omitnan'), 'Parent', a2);    
end

sb = repmat(ceil(sqrt(nReversals)), 1, 2);
savename1 = 'alwaysCsPlus_ch1_tiled';
h1 = ensureFigure(savename1, 1);
savename2 = 'alwaysCsPlus_ch2_tiled';
h2 = ensureFigure(savename2, 1);
for counter = 1:nReversals
    a1 = subplot(sb(1),sb(2),counter, 'Parent', h1); 
    plot(alwaysCsPlus_trialNumber, smoothdata(alwaysCsPlus_ch1(counter, :), 'movmean', 1, 'omitnan'), 'Parent', a1);
    a2 = subplot(sb(1),sb(2),counter, 'Parent', h2); 
    plot(alwaysCsPlus_trialNumber, smoothdata(alwaysCsPlus_ch2(counter, :), 'movmean', 1, 'omitnan'), 'Parent', a2);    
end

newCsPlus_licks = [AR.csMinus.csLicks.before AR.csPlus.csLicks.after];
newCsMinus_licks = [AR.csPlus.csLicks.before AR.csMinus.csLicks.after];
alwaysCsPlus_licks = [AR.csPlus.csLicks.before AR.csPlus.csLicks.after];


smoothWindow = 1;
newCsPlus_ch1_norm = smoothdata(newCsPlus_ch1, 2, 'movmean', smoothWindow, 'omitnan');
newCsPlus_ch2_norm = smoothdata(newCsPlus_ch2, 2, 'movmean', smoothWindow, 'omitnan');
newCsPlus_licks_norm = smoothdata(newCsPlus_licks, 2, 'movmean', smoothWindow, 'omitnan');

newCsMinus_ch1_norm = smoothdata(newCsMinus_ch1, 2, 'movmean', smoothWindow, 'omitnan');
newCsMinus_ch2_norm = smoothdata(newCsMinus_ch2, 2, 'movmean', smoothWindow, 'omitnan');
newCsMinus_licks_norm = smoothdata(newCsMinus_licks, 2, 'movmean', smoothWindow, 'omitnan');

alwaysCsPlus_ch1_norm = smoothdata(alwaysCsPlus_ch1, 2, 'movmean', smoothWindow, 'omitnan');
alwaysCsPlus_ch2_norm = smoothdata(alwaysCsPlus_ch2, 2, 'movmean', smoothWindow, 'omitnan');
alwaysCsPlus_licks_norm = smoothdata(alwaysCsPlus_licks, 2, 'movmean', smoothWindow, 'omitnan');

%% newCsPlus - normalize by f% of pre reversal values, smooth, find first and last common points across reversals

f = 0.5;
trialsBack = 30;
normVector_ch1 = percentile(newCsMinus_ch1_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), f, 2);
normVector_ch2 = percentile(newCsMinus_ch2_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), f, 2);
normVector_licks = percentile(newCsMinus_licks_norm(:,newCsPlus_firstRevTrial - trialsBack - 1:newCsPlus_firstRevTrial - 1), f, 2);
% normVector_ch1 = ones(nReversals, 1);
% normVector_ch2 = ones(nReversals, 1);
% normVector_licks = ones(nReversals, 1);

newCsPlus_ch1_norm = bsxfun(@rdivide, newCsPlus_ch1_norm, normVector_ch1);
newCsPlus_ch2_norm = bsxfun(@rdivide, newCsPlus_ch2_norm, normVector_ch2);
newCsPlus_licks_norm = bsxfun(@rdivide, newCsPlus_licks_norm, normVector_licks);



common = sum(~isnan(newCsPlus_ch1_norm)) > 3;


savename = 'reversals_newCsPlus';
ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_ch2_norm(:, common)), nanSEM(newCsPlus_ch2_norm(:, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_ch1_norm(:, common)), nanSEM(newCsPlus_ch1_norm(:, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(newCsPlus_licks_norm(:, common)), nanSEM(newCsPlus_licks_norm(:, common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
set(hla, 'LineWidth', 2);
% set(gca, 'XLim', [-40 40]);%, 'YLim', [-1 1]);
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



newCsMinus_ch1_norm = bsxfun(@rdivide, newCsMinus_ch1_norm, normVector_ch1);
newCsMinus_ch2_norm = bsxfun(@rdivide, newCsMinus_ch2_norm, normVector_ch2);
newCsMinus_licks_norm = bsxfun(@rdivide, newCsMinus_licks_norm, normVector_licks);


% common = find(mean(newCsMinus_ch1_norm));% applying mean will reveal NaNs
newCsMinus_common = sum(~isnan(newCsMinus_ch1_norm)) > 3;% applying mean will reveal NaNs

savename = 'reversals_newCsMinus';
ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmean(newCsMinus_ch2_norm(:, newCsMinus_common)), nanSEM(newCsMinus_ch2_norm(:, newCsMinus_common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmean(newCsMinus_ch1_norm(:, newCsMinus_common)), nanSEM(newCsMinus_ch1_norm(:, newCsMinus_common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(newCsMinus_common), nanmean(newCsMinus_licks_norm(:, newCsMinus_common)), nanSEM(newCsMinus_licks_norm(:, newCsMinus_common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
set(hla, 'LineWidth', 2);
% set(gca, 'XLim', [-40 40]);%, 'YLim', [-1 1]);
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

% alwaysCsPlus - normalize by f% of post reversal values, smooth, find first and last common points across reversals


alwaysCsPlus_ch1_norm = bsxfun(@rdivide, alwaysCsPlus_ch1_norm, normVector_ch1);
alwaysCsPlus_ch2_norm = bsxfun(@rdivide, alwaysCsPlus_ch2_norm, normVector_ch2);
alwaysCsPlus_licks_norm = bsxfun(@rdivide, alwaysCsPlus_licks_norm, normVector_licks);


% common = find(mean(alwaysCsPlus_ch1_norm));% applying mean will reveal NaNs
common = sum(~isnan(alwaysCsPlus_ch1_norm)) > 3;% applying mean will reveal NaNs


savename = 'reversals_alwaysCsPlus';
ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus_ch2_norm(:, common)), nanSEM(alwaysCsPlus_ch2_norm(:, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus_ch1_norm(:, common)), nanSEM(alwaysCsPlus_ch1_norm(:, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus_licks_norm(:, common)), nanSEM(alwaysCsPlus_licks_norm(:, common))',...
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