% LNL_odor_v2_pav_rev  reversal analysis script 9/
%% desktop 
files = {...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_17\', 'RE_DC_17.mat';...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_20\', 'RE_DC_20.mat';...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_35\', 'RE_DC_35.mat';...
    'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_36\', 'RE_DC_36.mat';...    
    };

%% laptop

files = {...
    'C:\Fitz_Data\SummaryAnalyses\LNL_Analysis_new\DC_17\', 'RE_DC_17.mat';...
    'C:\Fitz_Data\SummaryAnalyses\LNL_Analysis_new\DC_20\', 'RE_DC_20.mat';...
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
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SFN_2017\Reversals';
saveOn = 1;


%% find total number of reverals, max post reversal trials, etc.
% first column- CS+, second column, CS-
trialsBefore = zeros(length(RE), 2);
trialsAfter = zeros(length(RE), 2);
nReversals = zeros(length(RE), 1);
for counter = 1:length(RE)
    trialsBefore(counter, 1) = length(RE(counter).csPlus.trialsBefore);
    trialsAfter(counter, 1) = length(RE(counter).csPlus.trialsAfter);
    trialsBefore(counter, 2) = length(RE(counter).csMinus.trialsBefore);
    trialsAfter(counter, 2) = length(RE(counter).csMinus.trialsAfter);    
    nReversals(counter) = size(RE(1).csPlus.trialType.before, 1);
end


allReversals_phPeakMean_cs_ch1 = NaN(size();
allReversals_phPeakMean_cs_ch2 = NaN(sum(dataSize(:,1)), 82);

allReversals_phPeakMean_cs_ch1(1:6,1:82) = RE(1).csPlus.phPeakMean_cs_ch1.after(:,1:82);
allReversals_phPeakMean_cs_ch1(7:end,1:82) = RE(2).csPlus.phPeakMean_cs_ch1.after(:,1:82);


allReversals_phPeakMean_cs_ch2(1:6,1:82) = RE(1).csPlus.phPeakMean_cs_ch2.after(:,1:82);
allReversals_phPeakMean_cs_ch2(7:end,1:82) = RE(2).csPlus.phPeakMean_cs_ch2.after(:,1:82);

ensureFigure('Reversals_Raw', 1);
subplot(2,1,1); plot(allReversals_phPeakMean_cs_ch1');
subplot(2,1,2); plot(allReversals_phPeakMean_cs_ch2');

% tiled plots for ch1 and ch2
h1 = ensureFigure('raw_ch1_tiled', 1);
ylim1 = [min(min(allReversals_phPeakMean_cs_ch1)) max(max(allReversals_phPeakMean_cs_ch1))];
ylim2 = [min(min(allReversals_phPeakMean_cs_ch2)) max(max(allReversals_phPeakMean_cs_ch2))];
for counter = 1:10
    subplot(3,4,counter, 'Parent', h1); plot(allReversals_phPeakMean_cs_ch1(counter, :)); set(gca, 'YLim', ylim1);
end
h2 = ensureFigure('raw_ch2_tiled', 1);
for counter = 1:10
    subplot(3,4,counter, 'Parent', h2); plot(allReversals_phPeakMean_cs_ch2(counter, :)); set(gca, 'YLim', ylim2);
end

%% just plot the means
ensureFigure('means_raw', 1);
plot(nanmean(allReversals_phPeakMean_cs_ch1), 'g'); hold on
plot(nanmean(allReversals_phPeakMean_cs_ch2), 'r');

%% normalize by 90% of post reversal values

normVal_ch1 = zeros(10,1);
normVal_ch2 = zeros(10,1);
for counter = 1:10
    lastTrial = find(~isnan(allReversals_phPeakMean_cs_ch1(counter,:)), 1, 'last');
    normVal_ch1(counter) = percentile(allReversals_phPeakMean_cs_ch1(counter, :), 0.9);
    normVal_ch2(counter) = percentile(allReversals_phPeakMean_cs_ch2(counter, :), 0.9);
end

allReversals_phPeakMeanNorm_cs_ch1 = bsxfun(@rdivide, allReversals_phPeakMean_cs_ch1, normVal_ch1);
allReversals_phPeakMeanNorm_cs_ch2 = bsxfun(@rdivide, allReversals_phPeakMean_cs_ch2, normVal_ch2);

% reversal 6 is funky for both channels
ensureFigure('trials_norm', 1);
subplot(2,1,1); plot(allReversals_phPeakMeanNorm_cs_ch1'); set(gca, 'YLim', [-1 1.5]);
subplot(2,1,2); plot(allReversals_phPeakMeanNorm_cs_ch2'); set(gca, 'YLim', [-1 1.5]);

ensureFigure('means_norm', 1);
plot(nanmean(allReversals_phPeakMeanNorm_cs_ch1([1:5 7:10], :)), 'g'); hold on
plot(nanmean(allReversals_phPeakMeanNorm_cs_ch2([1:5 7:10], :)), 'g');
    
    
    
%% Same thing but add on before trials

%% 
% find total number of reverals, max post reversal trials, etc.
% hard coded for now
dataSize = [6 106; 4 82];
dataSizeBefore = [6 106; 4 80];

allReversals_phPeakMean_cs_ch1_before = NaN(sum(dataSizeBefore(:,1)), 80);
allReversals_phPeakMean_cs_ch2_before = NaN(sum(dataSizeBefore(:,1)), 80);

allReversals_phPeakMean_cs_ch1_before(1:6,1:80) = RE(1).csPlus.phPeakMean_cs_ch1.before(:,end - 80 + 1:end);
allReversals_phPeakMean_cs_ch1_before(7:end,1:80) = RE(2).csPlus.phPeakMean_cs_ch1.before(:,end - 80 + 1:end);


allReversals_phPeakMean_cs_ch2_before(1:6,1:80) = RE(1).csPlus.phPeakMean_cs_ch2.before(:,end - 80 + 1:end);
allReversals_phPeakMean_cs_ch2_before(7:end,1:80) = RE(2).csPlus.phPeakMean_cs_ch2.before(:,end - 80 + 1:end);

revNorm_ch1_before = bsxfun(@rdivide, allReversals_phPeakMean_cs_ch1_before, normVal_ch1);
revNorm_ch2_before = bsxfun(@rdivide, allReversals_phPeakMean_cs_ch2_before, normVal_ch2);

revNormComplete_ch1 = [revNorm_ch1_before allReversals_phPeakMeanNorm_cs_ch1];
revNormComplete_ch2 = [revNorm_ch2_before allReversals_phPeakMeanNorm_cs_ch2];


% tiled plots for ch1 and ch2
h1 = ensureFigure('normComplete_ch1_tiled', 1);
ylim1 = [min(min(revNormComplete_ch1)) max(max(revNormComplete_ch1))];
ylim2 = [min(min(revNormComplete_ch2)) max(max(revNormComplete_ch2))];

xData = -80:81;
for counter = 1:10
    subplot(3,4,counter, 'Parent', h1); plot(xData, revNormComplete_ch1(counter, end - 161:end)); set(gca, 'YLim', ylim1);
end

h2 = ensureFigure('normComplete_ch2_tiled', 1);
for counter = 1:10
    subplot(3,4,counter, 'Parent', h2); plot(xData, revNormComplete_ch2(counter, end - 161:end)); set(gca, 'YLim', ylim2);
end

%% by CS+
revNormComplete_ch1 = smoothdata(revNormComplete_ch1, 2, 'movmean', 3, 'omitnan');
revNormComplete_ch2 = smoothdata(revNormComplete_ch2, 2, 'movmean', 5, 'omitnan');
%%
% ensureFigure('full_means_norm', 1);
% plot(xData, nanmean(revNormComplete_ch1([1:5 7:10], :))), 'g'); hold on
% plot(xData, smooth(nanmean(revNormComplete_ch2([1:5 7:10], :))), 'r');
% set(gca, 'XLim', [-30 40]);
% {0.9258, 0.4883, 0.1914} % orange
% {0.6680, 0.2148, 0.8359} % purple
ensureFigure('Reversals_full_means_norm', 1);
hla = zeros(1,2);
[hl, hp] = boundedline(xData(18:127), nanmean(revNormComplete_ch2([1:5 7:10], 18:127)), nanSEM(revNormComplete_ch2([1:5 7:10], 18:127))',...
    'cmap', [237 125 49]/256, 'alpha'); hold on
hla(1) = hl;
[hl, hp] = boundedline(xData(18:127), nanmean(revNormComplete_ch1([1:5 7:10], 18:127)), nanSEM(revNormComplete_ch1([1:5 7:10], 18:127))',...
    'cmap', [171 55 214]/256, 'alpha');
hla(2) = hl;
set(hla, 'LineWidth', 2);
set(gca, 'XLim', [-20 40], 'YLim', [-1 1.2]);
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.'}, 'Location', 'southeast', 'FontSize', 18, 'Interpreter', 'tex', 'Box', 'off');

    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');
    
    formatFigureTalk([4 3]);
    if saveOn    
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm.fig'));
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm.jpg'));    
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm.emf'));                   
    end

ensureFigure('Reversals_full_means_norm_ChAT_only', 1);

[hl, hp] = boundedline(xData(18:127), nanmean(revNormComplete_ch1([1:5 7:10], 18:127)), nanSEM(revNormComplete_ch1([1:5 7:10], 18:127))',...
    'cmap', [171 55 214]/256, 'alpha');
set(hl, 'LineWidth', 2);
set(gca, 'XLim', [-20 40], 'YLim', [-1 1.2]);
    h  = addOrginLines;
    set(h, 'LineWidth', 2);
    legend(hl, {'\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.'}, 'Location', 'southeast', 'FontSize', 18, 'Interpreter', 'tex', 'Box', 'off');

    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');
    
    formatFigureTalk([4 3]);
    if saveOn    
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_ChAT_only.fig'));
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_ChAT_only.jpg'));    
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_ChAT_only.emf'));                   
    end


%% now track the new CS+ odor across the reversal

dataSizeBefore = [6 113; 4 73];

allReversals_phPeakMean_cs_ch1_before_byOdor = NaN(sum(dataSizeBefore(:,1)), 73);
allReversals_phPeakMean_cs_ch2_before_byOdor = NaN(sum(dataSizeBefore(:,1)), 73);

allReversals_phPeakMean_cs_ch1_before_byOdor(1:6,1:73) = RE(1).csMinus.phPeakMean_cs_ch1.before(:,end - 73 + 1:end);
allReversals_phPeakMean_cs_ch1_before_byOdor(7:end,1:73) = RE(2).csMinus.phPeakMean_cs_ch1.before(:,end - 73 + 1:end);


allReversals_phPeakMean_cs_ch2_before_byOdor(1:6,1:73) = RE(1).csMinus.phPeakMean_cs_ch2.before(:,end - 73 + 1:end);
allReversals_phPeakMean_cs_ch2_before_byOdor(7:end,1:73) = RE(2).csMinus.phPeakMean_cs_ch2.before(:,end - 73 + 1:end);

revNorm_ch1_before = bsxfun(@rdivide, allReversals_phPeakMean_cs_ch1_before_byOdor, normVal_ch1);
revNorm_ch2_before = bsxfun(@rdivide, allReversals_phPeakMean_cs_ch2_before_byOdor, normVal_ch2);

revNormComplete_ch1 = [revNorm_ch1_before allReversals_phPeakMeanNorm_cs_ch1];
revNormComplete_ch2 = [revNorm_ch2_before allReversals_phPeakMeanNorm_cs_ch2];


% tiled plots for ch1 and ch2
h1 = ensureFigure('normComplete_ch1_tiled', 1);
ylim1 = [min(min(revNormComplete_ch1)) max(max(revNormComplete_ch1))];
ylim2 = [min(min(revNormComplete_ch2)) max(max(revNormComplete_ch2))];

xData = -73:81;
for counter = 1:10
    subplot(3,4,counter, 'Parent', h1); plot(xData, revNormComplete_ch1(counter, :)); set(gca, 'YLim', ylim1);
end

h2 = ensureFigure('normComplete_ch2_tiled', 1);
for counter = 1:10
    subplot(3,4,counter, 'Parent', h2); plot(xData, revNormComplete_ch2(counter, :)); set(gca, 'YLim', ylim2);
end

% just plot the means
ensureFigure('full_means_norm_unrestricted', 1);
plot(xData, nanmean(revNormComplete_ch1([1:5 7:10], :)), 'g'); hold on
plot(xData, nanmean(revNormComplete_ch2([1:5 7:10], :)), 'r');


revNormComplete_delta = revNormComplete_ch1 - revNormComplete_ch2;
revNormComplete_ch1 = smoothdata(revNormComplete_ch1, 2, 'movmean', 3);
revNormComplete_ch2 = smoothdata(revNormComplete_ch2, 2, 'movmean', 5);
revNormComplte_delta = smoothdata(revNormComplete_delta, 2, 'movmean', 5);

%%
% {0.9258, 0.4883, 0.1914} % orange
% {0.6680, 0.2148, 0.8359} % purple
ensureFigure('Reversals_full_means_norm_ch1_only_byOdor', 1);
[hl, hp] = boundedline(xData(1:136), nanmean(revNormComplete_ch1([1:5 7:10], 1:136)), nanSEM(revNormComplete_ch1([1:5 7:10], 1:136))',...
    'alpha', 'cmap', [0.6680, 0.2148, 0.8359]); hold on
set(hl, 'LineWidth', 2);
% legend(hl, {'\color{green}ACh.'}, 'Location', 'southeast', 'FontSize', 18, 'Interpreter', 'tex', 'Box', 'off');
    set(gca, 'XLim', [-20 40], 'YLim', [0 1]);
    h  = addOrginLines;
    legend(hl, {'\bf\color[rgb]{0.6680, 0.2148, 0.8359}ACh.'}, 'Location', 'southeast', 'FontSize', 16, 'Interpreter', 'tex', 'Box', 'off');
    set(h, 'LineWidth', 2);
    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');
    
    formatFigureTalk([4 3]);
    if saveOn    
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_ch1_only_byOdor.fig'));
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_ch1_only_byOdor.jpg'));    
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_ch1_only_byOdor.meta'));                   
    end
ensureFigure('Reversals_full_means_norm_ch1_only_byOdor_matchedYaxis', 1);
[hl, hp] = boundedline(xData(1:136), nanmean(revNormComplete_ch1([1:5 7:10], 1:136)), nanSEM(revNormComplete_ch1([1:5 7:10], 1:136))',...
    'alpha', 'cmap', [0.6680, 0.2148, 0.8359]); hold on
set(hl, 'LineWidth', 2);
% legend(hl, {'\color{green}ACh.'}, 'Location', 'southeast', 'FontSize', 18, 'Interpreter', 'tex', 'Box', 'off');
    set(gca, 'XLim', [-20 40], 'YLim', [-1 1]);
    h  = addOrginLines;
    legend(hl, {'\bf\color[rgb]{0.6680, 0.2148, 0.8359}ACh.'}, 'Location', 'southeast', 'FontSize', 16, 'Interpreter', 'tex', 'Box', 'off');
    set(h, 'LineWidth', 2);
    xlabel('Odor presentations from reversal');
    ylabel('Cue response (norm.)');
    
    formatFigureTalk([4 3]);
    if saveOn    
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_ch1_only_byOdor_matchedYaxis.fig'));
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_ch1_only_byOdor_matchedYaxis.jpg'));    
        saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_ch1_only_byOdor_matchedYaxis.meta'));                   
    end

ensureFigure('Reversals_full_means_norm_byOdor', 1);

[lines, patches] = boundedline(xData(1:136), nanmean(revNormComplete_ch1([1:5 7:10], 1:136)), nanSEM(revNormComplete_ch1([1:5 7:10], 1:136))',...
    'alpha', 'cmap', [0.6680, 0.2148, 0.8359]); hold on
[hl, hp] = boundedline(xData(1:136), nanmean(revNormComplete_ch2([1:5 7:10], 1:136)), nanSEM(revNormComplete_ch2([1:5 7:10], 1:136))',...
    'cmap', [0.9258, 0.4883, 0.1914], 'alpha');
lines = [lines; hl];
set(lines, 'LineWidth', 2);
set(gca, 'XLim', [-20 40], 'YLim', [-1 1]);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(lines, {'\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.'}, 'Location', 'southeast', 'FontSize', 18, 'Interpreter', 'tex', 'Box', 'off');
xlabel('Odor presentations from reversal');
ylabel('Cue response (norm.)');

formatFigureTalk([4 3]);
if saveOn    
    saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_byOdor.fig'));
    saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_byOdor.jpg'));    
    saveas(gcf, fullfile(savepath, 'Reversals_full_means_norm_byOdor.meta'));                   
end
    
    
ensureFigure('full_means_normDelta', 1);
[hl, hp] = boundedline(xData(1:136), nanmean(revNormComplete_delta([1:5 7:10], 1:136)), nanSEM(revNormComplete_delta([1:5 7:10], 1:136))', 'k'); hold on
    set(gca, 'XLim', [-10 40]);%, 'YLim', [0 2]);
    














