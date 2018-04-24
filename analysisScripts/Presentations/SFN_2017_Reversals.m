
% SFN 2017 script 
% reversal data

%% desktop
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SFN_2017\Reversals';
saveOn = 1;

%%
figSize = [5 4];
fontSize = 24;
load('Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_20\TE.mat');
LNL_conditions;
% sorting by number of cue reward pairings post reversal

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
    [ha, hl, hp] = phPlotAverageFromTE(TE, trialMatrix(:,counter,:), 1,...
    'FluorDataField', 'ZS', 'window', [-1, 5], 'cmap', cmap(round(size(cmap, 1) * counter / nTranches), :), 'alpha', 0); %high value, reward
    patches(counter) = hp; lines(counter) = hl;
end
delete(patches); % WTF somehow not saving as vector file even though I have alpha set to off and I've tracked it all the way into boundedline
set(lines, 'LineWidth', 2);
ylabel('Fluor. (\sigma-baseline)'); set(gca, 'XLim', [-1, 5], 'YLim', [-.5 3]);
xlabel('Time from odor (s)');
% set(gca, 'XTick', []);
addStimulusPatch(gca, [0 1], '', patchHue, 1);
addStimulusPatch(gca, [2.9 3.1], '', patchHue, 1);
legend(lines, {'early', 'middle', 'late'}, 'Location', 'northwest', 'FontSize', 14, 'Box', 'off'); 
formatFigurePoster(figSize, '', fontSize);

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'epsc');
    saveas(gcf, fullfile(savepath, saveName), 'fig');   
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');             
end
    

% DAT
patchHue = [0.6 0.6 0.6];
saveName = 'rev_phAveragesByPairings_DAT';
ensureFigure(saveName, 1);
ax = axes; hold on; colormap(ax, 'autumn');
cmap = colormap(ax); % get the colormap
cmap = flipud(cmap);
cmap = cmap(1:round(size(cmap,1) * 0.66), :);
[ax lines patches] = deal(zeros(nTranches,1));
for counter = 1:nTranches
    [ha, hl, hp] = phPlotAverageFromTE(TE, trialMatrix(:,counter,:), 2,...
    'FluorDataField', 'ZS', 'window', [-1, 5], 'cmap', cmap(round(size(cmap, 1) * counter / nTranches), :), 'alpha', 0); %high value, reward
    patches(counter) = hp; lines(counter) = hl;
end
delete(patches); % WTF somehow not saving as vector file even though I have alpha set to off and I've tracked it all the way into boundedline
set(lines, 'LineWidth', 2);

ylabel('Fluor. (\sigma-baseline)'); set(gca, 'XLim', [-1, 5], 'YLim', [-.5 3]);
xlabel('Time from odor (s)');
addStimulusPatch(gca, [0 1], '', patchHue, 1);
addStimulusPatch(gca, [2.9 3.1], '', patchHue, 1);
legend(lines, {'early', 'middle', 'late'}, 'Location', 'northwest', 'FontSize', 14, 'Box', 'off'); 
formatFigurePoster(figSize, '', fontSize);

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'epsc');
    saveas(gcf, fullfile(savepath, saveName), 'fig');   
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');       
end
    

