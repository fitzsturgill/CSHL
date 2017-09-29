%% sorting by behavior (doesn't really work well)


excludeFirstTrials = 20;
excludeFirstTrials = TE.trialNumber > excludeFirstTrials;
threshPercentile = 0.15;

TE.csLicksAll = countEventFromTE(TE, 'Port1In', [0 3], TE.Cue);
ensureFigure('cueLick_Hist', 1);
[sorted index] = cum(TE.csLicksAll.rate(csPlusTrials & excludeFirstTrials));
plot(sorted, index); xlabel('lick rate');
threshRate = percentile(TE.csLicksAll.rate(csPlusTrials & excludeFirstTrials), threshPercentile);
hitTrialsRewarded = csPlusTrials & excludeFirstTrials & TE.csLicksAll.rate > threshRate & rewardTrials;
missTrialsRewarded = csPlusTrials & excludeFirstTrials & TE.csLicksAll.rate <= threshRate & rewardTrials;

ensureFigure('rev_phAveragesByBehavior', 1);
subplot(1, 2, 1);
[ha, hl] = phPlotAverageFromTE(TE, {hitTrialsRewarded, missTrialsRewarded}, 1,...
'FluorDataField', 'ZS', 'window', [-1, 5], 'linespec', {'g', 'r'}); %high value, reward
legend(hl, {'hit', 'miss'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
ylabel('CS+, outcomes'); set(gca, 'XLim', [-1, 5]); title('ACh.');%set(gca, 'YLim', ylim);

subplot(1,2,2);
[ha, hl] = phPlotAverageFromTE(TE, {hitTrialsRewarded, missTrialsRewarded}, 2,...
'FluorDataField', 'ZS', 'window', [-1, 5], 'linespec', {'g', 'r'}); %high value, reward
legend(hl, {'hit', 'miss'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
set(gca, 'XLim', [-1, 5]); title('Dop.')%set(gca, 'YLim', ylim);
formatFigureTalk([6,3]);

%% sorting by number of cue reward pairings post reversal

excludeFirstTrials = 20;
excludeFirstTrials = TE.trialNumber > excludeFirstTrials;
trancheSize = 7; % 
nSessions = max(TE.sessionIndex);
trialMatrix = zeros(length(TE.filename), nSessions);
trialMatrix2 = trialMatrix;
% find min number of reward/cue pairings across all reversals
minPairings = Inf;
for counter = 1:nSessions
    trialsThisSession = TE.sessionIndex == counter;
    thisReverse = find(TE.BlockChange & trialsThisSession);
    postReversalTrials = trialsThisSession; postReversalTrials(1:thisReverse - 1) = 0;
    pairingsThisReverse = postReversalTrials & csPlusTrials & rewardTrials;
    minPairings = min(minPairings, sum(pairingsThisReverse));
end
nTranches = floor(minPairings/trancheSize);
% subdivide cue/reward pairings trials into tranches for each
% session/reversal using reshape
trialMatrix = zeros(trancheSize, nTranches, nSessions);
for counter = 1:nSessions
    trialsThisSession = TE.sessionIndex == counter;
    thisReverse = find(TE.BlockChange & trialsThisSession);
    postReversalTrials = trialsThisSession; postReversalTrials(1:thisReverse - 1) = 0;
    pairingsThisReverse = postReversalTrials & csPlusTrials & rewardTrials;    
    pairingsThisReverse = find(pairingsThisReverse);
    trialMatrix(:,:,counter) = reshape(pairingsThisReverse(1:nTranches * trancheSize), [size(trialMatrix, 1), size(trialMatrix, 2)]);
end
% {0.9258, 0.4883, 0.1914} % orange
% {0.6680, 0.2148, 0.8359} % purple
%ChAT
patchHue = [0.6 0.6 0.6];
saveName = 'rev_phAveragesByPairings_ChAT';
ensureFigure(saveName, 1);
ax = axes; hold on; colormap(ax, 'cool');
cmap = colormap(ax); % get the colormap
cmap = cmap(round(size(cmap,1) * 0.33):end, :);
% cmap = flipud(cmap);
[patches lines] = deal(zeros(nTranches,1));
for counter = 1:nTranches
    [hp, hl] = phPlotAverageFromTE(TE, trialMatrix(:,counter,:), 1,...
    'FluorDataField', 'ZS', 'window', [-1, 5], 'cmap', cmap(round(size(cmap, 1) * counter / nTranches), :)); %high value, reward
    patches(counter) = hp; lines(counter) = hl;
end
set(lines, 'LineWidth', 2);
legend(lines, {'early', 'middle', 'late'}, 'Location', 'northwest', 'FontSize', 16, 'Box', 'off'); 
ylabel('\bf\color[rgb]{0.6680,0.2148,0.8359}Cholinergic \color{black}(\fontsize{20}\sigma\fontsize{16}-baseline)'); set(gca, 'XLim', [-1, 5], 'YLim', [-.5 3]);
% xlabel('Time from odor (s)');
set(gca, 'XTick', []);
addStimulusPatch(gca, [0 1], '', patchHue);
addStimulusPatch(gca, [2.9 3.1], '', patchHue);
formatFigureTalk;

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'epsc');
    saveas(gcf, fullfile(savepath, saveName), 'fig');   
    saveas(gcf, fullfile(savepath, saveName), 'meta');             
end
    

% DAT
patchHue = [0.6 0.6 0.6];
saveName = 'rev_phAveragesByPairings_DAT';
ensureFigure(saveName, 1);
ax = axes; hold on; colormap(ax, 'autumn');
cmap = colormap(ax); % get the colormap
cmap = flipud(cmap);
cmap = cmap(1:round(size(cmap,1) * 0.66), :);
[patches lines] = deal(zeros(nTranches,1));
for counter = 1:nTranches
    [hp hl] = phPlotAverageFromTE(TE, trialMatrix(:,counter,:), 2,...
    'FluorDataField', 'ZS', 'window', [-1, 5], 'cmap', cmap(round(size(cmap, 1) * counter / nTranches), :)); %high value, reward
    patches(counter) = hp; lines(counter) = hl;
end
set(lines, 'LineWidth', 2);
legend(lines, {'early', 'middle', 'late'}, 'Location', 'northwest', 'FontSize', 16, 'Box', 'off'); 
ylabel('\bf\color[rgb]{0.9258,0.4883,0.1914}Dopaminergic \color{black}(\fontsize{20}\sigma\fontsize{16}-baseline)'); set(gca, 'XLim', [-1, 5], 'YLim', [-.5 3]);
xlabel('Time from odor (s)');
addStimulusPatch(gca, [0 1], '', patchHue);
addStimulusPatch(gca, [2.9 3.1], '', patchHue);
formatFigureTalk;

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'epsc');
    saveas(gcf, fullfile(savepath, saveName), 'fig');   
    saveas(gcf, fullfile(savepath, saveName), 'meta');       
end
    

%% as above, outcome only

patchHue = [0.6 0.6 0.6];
saveName = 'rev_phAveragesByPairings_ChAT_outcomeOnly';
ensureFigure(saveName, 1);
window = [2 5];
figsize = [2.4 2.3];
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & periReversalTrials, csPlusTrials & rewardTrials & extraReversalTrials}, 1,...
'FluorDataField', 'ZS', 'window', window, 'linespec', {'b', 'k'}); %high value, reward
% legend(hl, {['\bf\color{cyan}'], '\bf\color{magenta}'}, 'Location', 'northoutside', 'FontSize', 16, 'Box', 'off'); 
ylabel('\bf\fontsize{20}\color{green} Cholinergic'); set(gca, 'XLim', window, 'YLim', [-.5 3]);
xlabel('Time from reward (s)');
addStimulusPatch(gca, [2.9 3.1], '', patchHue);
set(gca, 'XTickLabel', {'-1', '0', '1', '2'});
formatFigureTalk(figsize);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
end
    
% just crop legend in DAT version in powerpoint
saveName = 'rev_phAveragesByPairings_DAT_outcomeOnly';
ensureFigure(saveName, 1);
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & periReversalTrials, csPlusTrials & rewardTrials & extraReversalTrials}, 2,...
'FluorDataField', 'ZS', 'window', window, 'linespec', {'b', 'k'}); %high value, reward
% legend(hl, {['\bf\color{cyan} Initial trials'], '\bf\color{magenta} Others'}, 'Location', 'northoutside', 'FontSize', 16, 'Box', 'off'); 
set(gca, 'XLim', window, 'YLim', [-.5 3]);
ylabel('\bf\fontsize{20}\color{red} Dopaminergic');
xlabel('Time from reward (s)');
addStimulusPatch(gca, [2.9 3.1], '', patchHue);
set(gca, 'XTickLabel', {'-1', '0', '1', '2'});
formatFigureTalk(figsize);

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.emf']));   
end






















%% lets find optimal windows for CS for ACh and Dopamine

csPlus_special = csPlusTrials & excludeFirstTrials & TE.csLicksAll.rate > threshRate;
csMinus_special = csMinusTrials & excludeFirstTrials & TE.csLicksAll.rate <= threshRate;

ensureFigure('find_cue_windows', 1);
subplot(2,1,1);
[ha, hl] = phPlotAverageFromTE(TE, {csPlus_special, csMinus_special}, 1,...
'FluorDataField', 'ZS', 'window', [-.5, 1.5], 'linespec', {'g', 'r'}); %high value, reward
set(gca, 'XMinorTick', 'on', 'XGrid', 'on', 'XMinorGrid', 'on', 'GridColor', [0.3 0.3 0.3], 'MinorGridColor', [0.3 0.3 0.3], 'TickDir', 'both');
title('ACh');

subplot(2,1,2);
[ha, hl] = phPlotAverageFromTE(TE, {csPlus_special, csMinus_special}, 2,...
'FluorDataField', 'ZS', 'window', [-.5, 1.5], 'linespec', {'g', 'r'}); %high value, reward
set(gca, 'XMinorTick', 'on', 'XGrid', 'on', 'XMinorGrid', 'on', 'GridColor', [0.3 0.3 0.3], 'MinorGridColor', [0.3 0.3 0.3], 'TickDir', 'both');
title('Dopamine');

% set window to be 0.25 - 1 second after cue in LNL_odor_v2_pav_rev_AS