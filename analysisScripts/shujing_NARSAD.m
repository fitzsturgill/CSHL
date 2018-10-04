this MATLAB snippet: Untitled 
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\Grants\NARSAD_Shujing';

% savepath = 'C:\Users\Fitz\Dropbox\KepecsLab\_Fitz\Grants\NARSAD_Shujing';
load(fullfile('Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\ChAT_60', 'TE.mat'));

%% generate trial lookups for different combinations of conditions
saveOn = 0;
% reject session 4
% TE.reject(ismember(TE.sessionIndex, [1 2 3 4 5])) = 0;
% see Pavlovian_reversals_blocks  blocks 2 and 3
  validTrials = filterTE(TE, 'reject', 0);
  Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
  Odor2Trials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);  
  
  rewardTrials = filterTE(TE, 'trialType', [1 5], 'reject', 0);
  hitTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
  crTrials = filterTE(TE, 'trialOutcome', 2, 'reject', 0);  
  faTrials = filterTE(TE, 'trialOutcome', 0, 'reject', 0);  
  missTrials = filterTE(TE, 'trialOutcome', -1, 'reject', 0);
  punishTrials = filterTE(TE, 'trialType', [3 6], 'reject', 0);  
  neutralTrials = filterTE(TE, 'trialType', [2 4], 'reject', 0);
  block2Trials = filterTE(TE, 'BlockNumber', 2, 'reject', 0);
  block3Trials = filterTE(TE, 'BlockNumber', 3, 'reject', 0);
  csPlusTrials = filterTE(TE, 'trialType', [1 2], 'reject', 0);
  csMinusTrials = filterTE(TE, 'trialType', [3 4], 'reject', 0);
  uncuedReward = filterTE(TE, 'trialType', 5, 'reject', 0);
  uncuedPunish = filterTE(TE, 'trialType', 6, 'reject', 0);
  trialTypes = 1:6;
  trialsByType = cell(size(trialTypes));
  for counter = 1:length(trialTypes)
    trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
  end

trialCount = [1:length(TE.filename)]';


%% photometry averages, zscored
%   ylim = [-2 8];
  saveName = ['ChAT_60_BLA_phAvgs_late']; 
  h=ensureFigure(saveName, 1); 

  pm = [2 2];
  
  % - 6 0 4
 axes; 
set(gca, 'XLim', [-2 6], 'YLim', [-0.5 2]);  
addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7], 1); 
addStimulusPatch(gca, [2.9 3.1], '', [0.7 0.7 0.7], 1);
[ha, hl] = phPlotAverageFromTE(TE, {crTrials & punishTrials, crTrials & neutralTrials}, 1,...
  'FluorDataField', 'ZS', 'window', [-2, 6], 'linespec', {'r', 'k'}, 'alpha', ''); %high value, reward
%     legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');

ylabel('Fluor. (\fontsize{20}\sigma\fontsize{12}-baseline)'); xlabel('Time from odor (s)');  
  
  
  formatFigure('aspect', [1.2 1], 'scaleFactor', 2);
  if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');  
    saveas(gcf, fullfile(savepath, saveName), 'epsc');
  end  
  
%% raster plots
CLimFactor = 2;
% CLim = {[-0.005 0.005], [-0.1 0.1]};
% CLim = {[-0.005 0.005], [-0.1 0.1]};
trialStart = TE.Photometry.xData(1);
% reversals = find(TE.BlockChange);

  saveName = ['ChAT_60_BLA_phRasters_late'];
  h=ensureFigure(saveName, 1);
%   mcPortraitFigSetup(h);
  

  rasterTrials = crTrials & punishTrials & (TE.sessionIndex == 5);
  rasterTrials(897) = 0; % this one just looks weird
  axes; phRasterFromTE(TE, rasterTrials, 1, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive',...
    'medFilter', 5, 'window', [-2 6]);
    set(gca, 'FontSize', 10)
%   line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines  
%   colormap('jet');
    formatFigure('aspect', [1.2 1], 'scaleFactor', 2);
  if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');  
    saveas(gcf, fullfile(savepath, saveName), 'epsc');
  end  
    
    
    
    
%% something different, for ChAT_60 check if mouse is blinking with punishment

  %%
  behaviorTrials = crTrials & neutralTrials;
  behaviorTrials = csPlusTrials;
  ensureFigure('all_behavior_CsPlus', 1);
  reversals = find(diff(TE.BlockNumber(behaviorTrials, :))) + 1;
  subplot(1,4,1);
  
  image(TE.Wheel.data.V(behaviorTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%   set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
  set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
  line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines  
  title('Velocity');
  subplot(1,4,2);
  image(TE.pupil.pupDiameterNorm(behaviorTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
  set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
  line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines  
  colormap('parula'); 
  title('Pupil Diameter');  
  subplot(1,4,3); phRasterFromTE(TE, behaviorTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
  line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines  
  title('ChAT'); xlabel('Time frome odor (s)');
  subplot(1,4,4); phRasterFromTE(TE, behaviorTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
  line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines  
  title('DAT');  
  ax = findobj(gcf, 'Type', 'Axes');
%   set(ax, 'YLim', [40 80]);


%% early
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\Grants\NARSAD_Shujing';
CLimFactor = 2;
saveName = 'CuedOutcome_ChAT_60_phRasters_early';
  h = ensureFigure(saveName, 1); 
%   mcPortraitFigSetup(h);
  axes; phRasterFromTE(TE, trialsByType{5}, 1, 'CLimFactor', CLimFactor, 'medFilter', 5, 'window', [-5 3]);
%   title([TE.filename{1}(1:7) ': hival, punish'], 'Interpreter', 'none'); 
    formatFigure('aspect', [1.2 1], 'scaleFactor', 2);
  if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');  
    saveas(gcf, fullfile(savepath, saveName), 'epsc');
  end  
  
%% late

% find crap
crapTrials = find(crTrials & punishTrials & (TE.sessionIndex == 5), 30);
crapTrials = crapTrials(end - 10:end);
ensureFigure('findcrap', 1); 
for counter = 1:length(crapTrials)
  subplot(4,4, counter);
  plot(TE.Photometry.data.ZS(crap(counter), :))
end


%% expected uncertainty experiment
 load('Z:\SummaryAnalyses\LickNoLick_odor_v2_expected_uncertainty_all\DC_28\TE.mat');

% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks  blocks 2 and 3
validTrials = filterTE(TE, 'reject', 0);
highTrials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
mediumTrials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);  
lowTrials = filterTE(TE, 'OdorValveIndex', 3, 'reject', 0);  
rewardTrials = filterTE(TE, 'trialType', [1 3 5 7], 'reject', 0);
neutralTrials = filterTE(TE, 'trialType', [2 4 6], 'reject', 0);
trialTypes = 1:7;
trialsByType = cell(size(trialTypes));
for counter = 1:length(trialTypes)
  trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
end
uncuedReward = trialsByType{7};
trialCount = (1:length(TE.filename))';

% extract peak trial dFF responses to cues and reinforcement and lick counts
% csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);

cueWindow = [0 1];
traceWindow = [1 2];
usWindow = [0 0.5];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'

TE.cueLicks = countEventFromTE(TE, 'Port1In', cueWindow, TE.Cue);
TE.traceLicks = countEventFromTE(TE, 'Port1In', traceWindow, TE.Cue);
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);




TE.phPeakMean_baseline = bpCalcPeak_dFF(TE.Photometry, 1, [1 4], [], 'method', 'mean', 'phField', 'ZS');
TE.phPeakMean_us = bpCalcPeak_dFF(TE.Photometry, 1, usWindow, usZeros, 'method', 'mean', 'phField', 'ZS');

TE.phPeakMean_cue = bpCalcPeak_dFF(TE.Photometry, 1, cueWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
TE.phPeakPercentile_cue = bpCalcPeak_dFF(TE.Photometry, 1, cueWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');

TE.phPeakMean_trace = bpCalcPeak_dFF(TE.Photometry, 1, traceWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
TE.phPeakPercentile_trace = bpCalcPeak_dFF(TE.Photometry, 1, traceWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');

TE.phPeakPercentile_us = bpCalcPeak_dFF(TE.Photometry, 1, [0 0.75], usZeros, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');

%%
saveOn = 1;
saveName = 'expUncertainty_phAverage';  
h=ensureFigure(saveName, 1); 
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\Grants\NARSAD_Shujing';
% savepath = 'C:\Users\Fitz\Dropbox\KepecsLab\_Fitz\Grants\NARSAD_Shujing';




% First row: by cue condition
% varargin = {'trialNumbering', 'consecutive',...
%   'window', [-3 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
% axh = [];
% subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, {highTrials, mediumTrials, lowTrials}, 'Port1In', varargin{:});
% legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
% title('Licking'); ylabel('Cue only'); textBox(TE.filename{1}(1:6)); set(gca, 'XLim', [-3 3]);

axes; 
set(gca, 'XLim', [-2 2], 'YLim', [- 1 3]);
addStimulusPatch(gca, [0 1], '', [0.9 0.9 0.9], 1); 
addStimulusPatch(gca, [1 2], '', [0.5 1 0.5], 1); 
[ha, hl] = phPlotAverageFromTE(TE, {highTrials, mediumTrials, lowTrials, uncuedReward}, 1, 'FluorDataField', 'ZS', 'window', [-2 2], 'alpha', [],...
  'linespec', {'c', 'b', 'm', 'k'}); %high value, reward
legend(hl, {'high', 'medium', 'low', 'uncued'}, 'Location', 'northwest', 'FontSize', 8); legend('boxoff');
xlabel('Time from odor (s)'); ylabel('Fluor. (\fontsize{20}\sigma\fontsize{12}-baseline)'); 
formatFigure('aspect', [1.2 1], 'scaleFactor', 2);
if saveOn
  saveas(gcf, fullfile(savepath, saveName), 'fig');
  saveas(gcf, fullfile(savepath, saveName), 'jpeg');  
  saveas(gcf, fullfile(savepath, saveName), 'epsc');
end  

%% 
trialTypes = {lowTrials; mediumTrials; highTrials};

chat.avg = zeros(length(trialTypes), 2);
chat.sem = zeros(length(trialTypes), 2);

% somewhat kludgy way of baselining and normalizing data prior to determing
% mean and SEM
dataAvgs = zeros(length(trialTypes), 2);
for counter = 1:length(trialTypes)
  dataAvgs(counter, 1) = mean(TE.phPeakMean_cue.data(trialTypes{counter}));
  dataAvgs(counter, 2) = mean(TE.phPeakMean_trace.data(trialTypes{counter}));
end
%

for counter = 1:length(trialTypes)
  data = TE.phPeakMean_cue.data(trialTypes{counter});
  data = data - dataAvgs(1, 1); % baseline
  data = data / dataAvgs(end, 1); % normalize
  chat.avg(counter, 1) = mean(data);
  chat.sem(counter, 1) = std(data) ./ sqrt(length(data));
  data = TE.phPeakMean_trace.data(trialTypes{counter});
  data = data - dataAvgs(1, 2);
  data = data / dataAvgs(end, 2);
  chat.avg(counter, 2) = mean(data);
  chat.sem(counter, 2) = std(data) ./ sqrt(length(data));  
end

chat.avg = bsxfun(@minus, chat.avg, chat.avg(1,:));
chat.avg = bsxfun(@rdivide, chat.avg, chat.avg(end,:));
saveName = 'uncertainty_summary';
h = ensureFigure(saveName, 1);

axes;
b = errorbar(chat.avg, chat.sem, 'LineWidth', 1);
b(1).Color = [0.5 0.5 0.5];
b(2).Color = [0 1 0];
% ylabel('\bf\color[rgb]{0.6680,0.2148,0.8359}Cholinergic \color{black}(\fontsize{20}\sigma\fontsize{16}-baseline)');
% legend(b, {'\color[rgb]{0.5 1 0.5} odor', '\color[rgb]{0 1 0} trace'}, 'Location', 'northwest', 'FontSize', 10, 'Interpreter', 'tex', 'Box', 'off');
legend(b, {'odor', 'trace'}, 'Location', 'northwest', 'FontSize', 10, 'Box', 'off');
% legend(b, {['\bf\color{cyan}'], '\bf\color{magenta}'}, 'Location', 'northoutside', 'FontSize', 16, 'Box', 'off'); 
set(gca, 'XLim', [0.5 3.5], 'XTick', [1 2 3], 'XTickLabel', {'0.25', '0.5', '0.75'});
xlabel('Probability of Reward'); ylabel('Fluor. (normalized)');

formatFigure('aspect', [1.2 1], 'scaleFactor', 2);
if saveOn
  saveas(gcf, fullfile(savepath, saveName), 'fig');
  saveas(gcf, fullfile(savepath, saveName), 'jpeg');  
  saveas(gcf, fullfile(savepath, saveName), 'epsc');
end  