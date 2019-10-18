% reversals_noPunish_exampleMouse

% using... DC_56



% value reversal
DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);
smoothWindow = 1;
saveOn = 1;




animals = {'DC_56', 'DC_35'};
% animals = {'DC_35'};
sessions = [1 3];
% sessions = [3];

photometryField = 'Photometry';
% photometryField = 'PhotometryExpFit';
fdField = 'ZS';


    
%% plot odor 1 and odor 2 Rasters for an example session
    figSize = [4 4];
% Define the positions of different axes matrices on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    params = struct();    
    params.matpos = [0 0 1 1];    
%     params.cellmargin = [.05 .05 0.05 0.05];    
    params.figmargin = [0.15 0 0.1 0.1];
    params.cellmargin = [0.025 0.025 0.025 0.025];    
    
for counter = 1:length(animals)
    animal = animals{counter};
    session = sessions(counter);
    success = dbLoadAnimal(DB, animal); % load TE and trial lookups

    saveName = sprintf('example_allBehavior_v2_%s', animal); 
    ensureFigure(saveName, 1);

    hax = axesmatrix(2,3,1:6,params);
    exampleTrials = Odor1Trials & ismember(TE.sessionIndex, session);
    reversals = find(diff(TE.BlockNumber(exampleTrials, :))) + 1;
    sessionChanges = find(diff(TE.sessionIndex(exampleTrials, :))) + 1;

    climfactor = 3;
    axes(hax(1));
    [~, lh] = eventRasterFromTE(TE, exampleTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7], 'XTick', []);
    set(lh, 'LineWidth', 0.3, 'Color', [0 0 0]);
    title('licking');
    ylabel('Odor 1 trial number');

    axes(hax(2));
    phRasterFromTE(TE, exampleTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor,...
        'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    tcolor = mycolors('chat');
    set(gca, 'XTick', [], 'YTick', []);
    title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Ach.']);
    
    axes(hax(3));
    phRasterFromTE(TE, exampleTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor,...
        'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    tcolor = mycolors('dat');
    set(gca, 'XTick', [], 'YTick', []);
    title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Dop.']);

    exampleTrials = Odor2Trials & ismember(TE.sessionIndex, session);
    reversals = find(diff(TE.BlockNumber(exampleTrials, :))) + 1;
    sessionChanges = find(diff(TE.sessionIndex(exampleTrials, :))) + 1;

    axes(hax(4));
    [~, lh] = eventRasterFromTE(TE, exampleTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7]);
    set(lh, 'LineWidth', 0.3, 'Color', [0 0 0]);
    ylabel('Odor 2 trial number');

    axes(hax(5));
    phRasterFromTE(TE, exampleTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor,...
        'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    xlabel('Time from odor (s)');
    set(gca, 'YTick', []);
    
    axes(hax(6));
    phRasterFromTE(TE, exampleTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor,...
        'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    set(gca, 'YTick', []);    




    formatFigurePublish('size', figSize);
    if saveOn 
        print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%         export_fig(fullfile(savepath, saveName), '-eps');
    end
end

%% averages

animal = 'DC_56';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups
fdField = 'ZS';
saveName = sprintf('%s_phAvgs_%s', animal, fdField);  
h=ensureFigure(saveName, 1); 


linecolors = [1 0 0; 0 0 1; 0 1 1; 0 1 0];            

subplot(1, 2, 1);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Ach.']);
addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [2.9 3.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {neutralTrials & csPlusTrials & hitTrials, rewardTrials & csPlusTrials & hitTrials, csMinusTrials & CRTrials & rewardTrials, uncuedReward}, 1,...
    'FluorDataField', fdField, 'window', [-4, 7], 'cmap', linecolors); %high value, reward

% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northwest'); legend('boxoff');
% ylabel('\fontsize{8}Fluor. (\fontsize{12}\sigma\fontsize{8}-bl.)'); 
xlabel('time from cue (s)'); set(gca, 'XLim', [-4 7]);



subplot(1, 2, 2);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('dat');
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Dop.']);
addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [2.9 3.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {neutralTrials & csPlusTrials & hitTrials, rewardTrials & csPlusTrials & hitTrials, csMinusTrials & CRTrials & rewardTrials, uncuedReward}, 2,...
    'FluorDataField', fdField, 'window', [-4, 7], 'cmap', linecolors); %high value, reward
% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northoutside'); legend('boxoff');
% ylabel('\fontsize{8}Fluor. (\fontsize{12}\sigma\fontsize{8}-bl.)'); 
xlabel('time from cue (s)'); set(gca, 'XLim', [-4 7]);

formatFigurePublish('size', [2 1]);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% SIGNAL CORRELATIONS: show that cue and reward responses are correlated on a trial-by-trial basis between CBF and VTA
linecolors = [1 0 0; 0 0 1; 0 1 1; 0 1 0];         
saveName = 'CBF_VTA_correlations_reversals';
ensureFigure(saveName, 1);
trialSets = [neutralTrials & csPlusTrials & hitTrials, rewardTrials & csPlusTrials & hitTrials, csMinusTrials & CRTrials & rewardTrials, uncuedReward];
allTrials = sum(trialSets, 2) ~= 0;
xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];
% ylims = [min(TE.phPeakMean_cs(1).data(allTrials) max(TE.phPeakMean_cs(1).data(allTrials); min(TE.phPeakMean_us(1).data(allTrials) max(TE.phPeakMean_us(1).data(allTrials)];

% allTrials = true; 
subplot(1,2,1); hold on; subplot(1,2,2); hold on;
for counter = 1:size(trialSets, 2)
    subplot(1,2,1); set(gca, 'XLim', xlims(1,:));
%     allTrials = allTrials | trialSets{counter};
    xData = TE.phPeakMean_cs(2).data(trialSets(:, counter)); yData = TE.phPeakMean_cs(1).data(trialSets(:, counter)); 
%     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for cs
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off; %,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));

    subplot(1,2,2); set(gca, 'XLim', xlims(2,:));
    xData = TE.phPeakMean_us(2).data(trialSets(:, counter)); yData = TE.phPeakMean_us(1).data(trialSets(:, counter)); 
%     scatter(TE.phPeakMean_us(2).data(trialSets{counter}), TE.phPeakMean_us(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');    
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for us
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off;% ,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
end
% % fit for cs
subplot(1,2,1);
% xData = TE.phPeakMean_cs(2).data(allTrials); yData = TE.phPeakMean_cs(1).data(allTrials); 
% fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
% fob = fit(xData, yData, 'poly1', fo); 
% fph=plot(fob,'predfunc'); legend off;
% set(fph, 'LineWidth', 0.5, 'Color', 'k');
% % textBox(['R = ' num2str(corr(xData, yData), 3)],[], [0.4 0.95], 8);
textBox('Cue', [], [0.5 0.95], 8);
xlabel('\fontsize{8}Dop. (\fontsize{12}\sigma\fontsize{8}-baseline)');
ylabel('\fontsize{8}Ach. (\fontsize{12}\sigma\fontsize{8}-baseline)');
% 
% % fit for us
subplot(1,2,2);
% xData = TE.phPeakMean_us(2).data(allTrials); yData = TE.phPeakMean_us(1).data(allTrials); 
% fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
% fob = fit(xData, yData, 'poly1', fo); 
% fph=plot(fob,'predfunc'); legend off;
% set(fph, 'LineWidth', 0.5, 'Color', 'k');
% % textBox(['R = ' num2str(corr(xData, yData), 3)],[], [0.4 0.95], 8);
textBox('Reward', [], [0.5 0.95], 8);
xlabel('\fontsize{8}Dop. (\fontsize{12}\sigma\fontsize{8}-baseline)');
ylabel('');

formatFigurePublish('size', [2.5 1.1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% NOISE CORRELATIONS: show that cue and reward responses are correlated on a trial-by-trial basis between CBF and VTA
linecolors = [1 0 0; 0 0 1; 0 1 1; 0 1 0];         
saveName = 'CBF_VTA_NOISEcorrelations_reversals';
ensureFigure(saveName, 1);
trialSets = [neutralTrials & csPlusTrials & hitTrials, rewardTrials & csPlusTrials & hitTrials, csMinusTrials & CRTrials & rewardTrials, uncuedReward];
allTrials = sum(trialSets, 2) ~= 0;
xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];
ylims = [min(TE.phPeakMean_cs(1).data(allTrials)) max(TE.phPeakMean_cs(1).data(allTrials)); min(TE.phPeakMean_us(1).data(allTrials)) max(TE.phPeakMean_us(1).data(allTrials))];

% allTrials = true; 
subplot(1,2,1); hold on; subplot(1,2,2); hold on;
for counter = 1:size(trialSets, 2)
    subplot(1,2,1); set(gca, 'XLim', xlims(1,:)); set(gca, 'YLim', ylims(1,:));
%     allTrials = allTrials | trialSets{counter};
    xData = TE.phPeakMean_cs(2).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_cs(2).data(trialSets(:, counter))); yData = TE.phPeakMean_cs(1).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_cs(1).data(trialSets(:, counter))); 
%     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for cs
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off; %,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));

    subplot(1,2,2); set(gca, 'XLim', xlims(2,:)); set(gca, 'YLim', ylims(2,:));
    xData = TE.phPeakMean_us(2).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_us(2).data(trialSets(:, counter))); yData = TE.phPeakMean_us(1).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_us(1).data(trialSets(:, counter))); 
%     scatter(TE.phPeakMean_us(2).data(trialSets{counter}), TE.phPeakMean_us(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');    
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for us
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off;% ,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
end
% % fit for cs
subplot(1,2,1);
% xData = TE.phPeakMean_cs(2).data(allTrials); yData = TE.phPeakMean_cs(1).data(allTrials); 
% fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
% fob = fit(xData, yData, 'poly1', fo); 
% fph=plot(fob,'predfunc'); legend off;
% set(fph, 'LineWidth', 0.5, 'Color', 'k');
% % textBox(['R = ' num2str(corr(xData, yData), 3)],[], [0.4 0.95], 8);
textBox('Cue', [], [0.5 0.95], 8);
xlabel('\fontsize{8}Dop. (\fontsize{12}\sigma\fontsize{8}-baseline)');
ylabel('\fontsize{8}Ach. (\fontsize{12}\sigma\fontsize{8}-baseline)');
% 
% % fit for us
subplot(1,2,2);
% xData = TE.phPeakMean_us(2).data(allTrials); yData = TE.phPeakMean_us(1).data(allTrials); 
% fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
% fob = fit(xData, yData, 'poly1', fo); 
% fph=plot(fob,'predfunc'); legend off;
% set(fph, 'LineWidth', 0.5, 'Color', 'k');
% % textBox(['R = ' num2str(corr(xData, yData), 3)],[], [0.4 0.95], 8);
textBox('Reward', [], [0.5 0.95], 8);
xlabel('\fontsize{8}Dop. (\fontsize{12}\sigma\fontsize{8}-baseline)');
ylabel('');

formatFigurePublish('size', [2.5 1.1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% NOISE CORRELATIONS #2: cue and reward combined on one figure
linecolors = [1 0 0; 0 0 1; 0 1 1; 0 1 0];       
markerSize = 4;
figSize = [1 1];
saveName = 'CBF_VTA_NOISEcorrelations_reversals_v2';
ensureFigure(saveName, 1);
trialSets = [neutralTrials & csPlusTrials & hitTrials, rewardTrials & csPlusTrials & hitTrials, csMinusTrials & CRTrials & rewardTrials, uncuedReward];
allTrials = sum(trialSets, 2) ~= 0;
xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];
ylims = [min(TE.phPeakMean_cs(1).data(allTrials)) max(TE.phPeakMean_cs(1).data(allTrials)); min(TE.phPeakMean_us(1).data(allTrials)) max(TE.phPeakMean_us(1).data(allTrials))];

% allTrials = true; 
axes; hold on;
set(gca, 'XLim', xlims(1,:)); set(gca, 'YLim', ylims(1,:));
for counter = 1:size(trialSets, 2)

%     allTrials = allTrials | trialSets{counter};
    xData = TE.phPeakMean_cs(2).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_cs(2).data(trialSets(:, counter))); yData = TE.phPeakMean_cs(1).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_cs(1).data(trialSets(:, counter))); 
%     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
    scatter(xData, yData, markerSize, linecolors(counter, :), '.');
    % fit for cs
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off; %,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :), 'LineStyle', '--');
    
    xData = TE.phPeakMean_us(2).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_us(2).data(trialSets(:, counter))); yData = TE.phPeakMean_us(1).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_us(1).data(trialSets(:, counter))); 
%     scatter(TE.phPeakMean_us(2).data(trialSets{counter}), TE.phPeakMean_us(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');    
    scatter(xData, yData, markerSize, linecolors(counter, :), '.');
    % fit for us
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off;% ,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
end
setXYsameLimit;
xlabel('\fontsize{8}Dop. (\fontsize{12}\sigma\fontsize{8}-bl)');
ylabel('\fontsize{8}Ach. (\fontsize{12}\sigma\fontsize{8}-bl)');
formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end






%% plot the different quintiles of csPlusReward trials conditioned on csPlus cue photometry quintile

split = percentile(TE.phPeakMean_cs(1).data(csPlusTrials), 0.25);
subplot(2, 2, 3);
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & (TE.phPeakMean_cs(1).data < split), csPlusTrials & rewardTrials & (TE.phPeakMean_cs(1).data >= split)}, 1,...
'FluorDataField', fdField, 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
title('CS+, outcomes'); set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);

split = percentile(TE.phPeakMean_cs(1).data(csPlusTrials), 0.25);

subplot(2, 2, 4);
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & (TE.phPeakMean_cs(2).data < split), csPlusTrials & rewardTrials & (TE.phPeakMean_cs(2).data >= split)}, 2,...
    'FluorDataField', fdField, 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
xlabel('time from cue (s)');     set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);

%% try conditioning on cue response in ch1 and 2



%% cumulative histogram of reward responses conditioned on behavior (but not cue condition


bottom = (TE.licks_cs.rate == 0);
top = (TE.licks_cs.rate > 0);
bottom_ch1 = cum(TE.phPeakPercentile_us(1).data(bottom & rewardTrials));
top_ch1 = cum(TE.phPeakPercentile_us(1).data(top & rewardTrials));
bottom_ch2 = cum(TE.phPeakPercentile_us(2).data(bottom & rewardTrials));
top_ch2 = cum(TE.phPeakPercentile_us(2).data(top & rewardTrials));

% first the averages
saveName = 'Reward_cumHist_companionAverages';
ensureFigure(saveName, 1); 

[~, hl, hp] = phPlotAverageFromTE(TE, {bottom & rewardTrials, top & rewardTrials}, 1,...
    'FluorDataField', fdField, 'window', [-1, 6], 'cmap', repmat(mycolors('chat'), 2, 1)); %high value, reward
set(hl(2), 'LineStyle', ':', 'LineWidth', 1);
set(hl(1), 'LineStyle', '-', 'LineWidth', 1);
set(hp(2), 'FaceAlpha', 0.1);
hold on;
[~, hl, hp] = phPlotAverageFromTE(TE, {bottom & rewardTrials, top & rewardTrials}, 2,...
    'FluorDataField', fdField, 'window', [-1, 6], 'cmap', repmat(mycolors('dat'), 2, 1)); %high value, reward
set(hl(2), 'LineStyle', ':', 'LineWidth', 1);
set(hl(1), 'LineStyle', '-', 'LineWidth', 1);
set(hp(2), 'FaceAlpha', 0.1);
xlabel('time from cue (s)'); ylabel('\fontsize{8}Fluor. (\fontsize{12}\sigma\fontsize{8}-baseline)');
yMin = get(gca, 'YLim'); yMin = yMin(1);
set(gca, 'YLim', [yMin 8]);    
addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [2.9 3.1], '', [0.7 0.7 0.7], 0.4);

formatFigurePublish('size', [1.2 1.2]);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end


saveName = 'Reward_cumHist';
ensureFigure(saveName, 1);
plot(bottom_ch1.sorted, bottom_ch1.index, '-', 'Color', mycolors('chat')); hold on;
plot(top_ch1.sorted, top_ch1.index, ':', 'Color', mycolors('chat')); 
plot(bottom_ch2.sorted, bottom_ch2.index, '-', 'Color', mycolors('dat')); 
plot(top_ch2.sorted, top_ch2.index, ':', 'Color', mycolors('dat'));
xlabel('\fontsize{8}Fluor. (\fontsize{12}\sigma\fontsize{8}-baseline (mean))');
% legend({'', '', '', ''}, 'Location', 'southeast'); legend('boxoff');
set(gca, 'XLim', [-5 20]);

formatFigurePublish('size', [1.2 1.2]);
if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

% subplot(1,2,2);
% % scatter(TE.licks_cs.count(rewardTrials), TE.phPeakMean_us(2).data(rewardTrials), '.');
% lickBins = [0:5 10] - 0.5;
% [ch1_means, ch1_sem, ~] = binnedMeansXY(TE.licks_cs.count(rewardTrials), TE.phPeakMean_us(1).data(rewardTrials), lickBins);
% [ch2_means, ch2_sem, binCenters] = binnedMeansXY(TE.licks_cs.count(rewardTrials), TE.phPeakMean_us(2).data(rewardTrials), lickBins);
% errorbar(binCenters, ch1_means, ch1_sem, 'Color', mycolors('chat')); hold on;
% errorbar(binCenters, ch2_means, ch2_sem, 'Color', mycolors('dat'));
% if saveOn
%     saveas(gcf, fullfile(savepath, [saveName '.fig']));
%     saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
% end    




%% Reversal Averages


% load the reversals
load(fullfile(DB.path, 'pooled', sprintf('RE_%s.mat', animal)));    
%% find good reversals
% best reversal is reversal #1
smoothWindow = 3;
% plot individual reversals
nRev = size(RE.csPlus.phPeakMean_cs_ch1.before, 1);
ensureFigure('findGoodRevs', 1);
xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
% bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);
% ceilIx(1) = nearest(RE.csMinus.trialsBefore, -30); ceilIx(2) = nearest(RE.csPlus.trialsBefore, 0);
for counter = 1:nRev
    subplot(3,4,counter); hold on;
    plot(xData, smoothdata([RE.csMinus.licks_cs.before(counter,:) RE.csPlus.licks_cs.after(counter, :)], 2, 'movmean', smoothWindow),'k'); hold on;
%     plot(xData, smoothdata([RE.csMinus.phPeakMean_cs_ch1.before(counter,:) RE.csPlus.phPeakMean_cs_ch1.after(counter, :)], 2, 'movmean', smoothWindow), 'g');
%     plot(xData, smoothdata([RE.csMinus.phPeakMean_cs_ch2.before(counter,:) RE.csPlus.phPeakMean_cs_ch2.after(counter, :)], 2, 'movmean', smoothWindow), 'r');
%     plot(xData, smoothdata([RE.csMinus.licks_cs.before(counter,:) RE.csPlus.licks_cs.after(counter, :)], 2, 'movmean', 1),'k.'); hold on;
    plot(xData, smoothdata([RE.csMinus.phPeakMean_cs_ch1.before(counter,:) RE.csPlus.phPeakMean_cs_ch1.after(counter, :)], 2, 'movmean', 1), 'g.');
    plot(xData, smoothdata([RE.csMinus.phPeakMean_cs_ch2.before(counter,:) RE.csPlus.phPeakMean_cs_ch2.after(counter, :)], 2, 'movmean', 1), 'r.');
    textBox(int2str(counter));
end
%% CS+: plot reversal # 1 and fit a weibull function to it
smoothWindow = 3;
% plot individual reversals
nRev = size(RE.csPlus.phPeakMean_cs_ch1.before, 1);
rev = 1;
saveName = 'exampleReversal_newCsPlus';
ensureFigure(saveName, 1);
ch1Data = [RE.csMinus.phPeakMean_cs_ch1.before(rev,:) RE.csPlus.phPeakMean_cs_ch1.after(rev, :)];
ch2Data = [RE.csMinus.phPeakMean_cs_ch2.before(rev,:) RE.csPlus.phPeakMean_cs_ch2.after(rev, :)];
lickData = [RE.csMinus.licks_cs.before(rev,:) RE.csPlus.licks_cs.after(rev, :)];
validtrials = ~isnan(ch1Data);
ch1Data = ch1Data(validtrials);
ch2Data = ch2Data(validtrials);
lickData = lickData(validtrials);

xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter] + 1; % doesn't make sense that last trial of last reversal is trial 0 does it?
xData = xData(validtrials);
model = 'a * (1 - exp(-1 * (x/b)^c)) + d'; % weibull function, CDF form
fo = fitoptions('Method', 'NonlinearLeastSquares',...
    'Upper', [Inf  Inf Inf Inf],...
    'Lower', [0 0 0 -Inf],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
    'StartPoint', [1 round(length(ch1Data)/2) 1 0]...
    );
ft = fittype(model, 'options', fo);
%     xData = (0:length(toFit) - 1)';
[fitobject_ch1, ~, ~] = fit((0:(length(ch1Data) - 1))', ch1Data', ft, fo);
[fitobject_ch2, ~, ~] = fit((0:(length(ch1Data) - 1))', ch2Data', ft, fo);

axes; yyaxis top; hold on;
hl = [];
hl(end + 1) = plot(xData, ch1Data, '.', 'Color', mycolors('chat'), 'MarkerSize', 4);
hl(end + 1) = plot(xData, ch2Data, '.', 'Color', mycolors('dat'), 'MarkerSize', 4);
plot(xData, fitobject_ch1.a .* (1 - exp(-1 * ((xData - xData(1)) ./ fitobject_ch1.b) .^fitobject_ch1.c)) + fitobject_ch1.d, '--', 'Color', mycolors('chat'));
plot(xData, fitobject_ch2.a .* (1 - exp(-1 * ((xData - xData(1)) ./ fitobject_ch2.b) .^fitobject_ch2.c)) + fitobject_ch2.d, '--', 'Color', mycolors('dat'));
set(gca, 'YLim', [-4 7]);    
yyaxis right; set(gca, 'YColor', [0 0 0]); ylabel('licks/s');
hl(end + 1) = plot(xData, smoothdata(lickData, 2, 'movmedian', 5),'k');
set(gca, 'YLim', [0 6], 'XLim', xData([1 end]));
yyaxis left;
xlabel('New Cs+ trials from reversal')
ylabel('\fontsize{8}Fluor. (\fontsize{12}\sigma\fontsize{8}-baseline)');
% legend({'ACh.', 'Dop.', 'Licks'}, 'Box', 'off', 'Location', 'best', 'Orientation', 'vertical');
formatFigurePublish('size', [1.5 1.2]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% CS-: plot reversal # 1 and fit a weibull function to it
smoothWindow = 3;
% plot individual reversals
nRev = size(RE.csMinus.phPeakMean_cs_ch1.before, 1);
rev = 1;
saveName = 'exampleReversal_newCsMinus';
ensureFigure(saveName, 1);
ch1Data = [RE.csPlus.phPeakMean_cs_ch1.before(rev,:) RE.csMinus.phPeakMean_cs_ch1.after(rev, :)];
ch2Data = [RE.csPlus.phPeakMean_cs_ch2.before(rev,:) RE.csMinus.phPeakMean_cs_ch2.after(rev, :)];
lickData = [RE.csPlus.licks_cs.before(rev,:) RE.csMinus.licks_cs.after(rev, :)];
validtrials = ~isnan(ch1Data);
ch1Data = ch1Data(validtrials);
ch2Data = ch2Data(validtrials);
lickData = lickData(validtrials);

xData = [RE.csPlus.trialsBefore RE.csMinus.trialsAfter] + 1; % doesn't make sense that last trial of last reversal is trial 0 does it?
xData = xData(validtrials);
model = 'a * (exp(-1 * (x/b)^c)) + d'; % negative slope version of weibull function, CDF form
fo = fitoptions('Method', 'NonlinearLeastSquares',...
    'Upper', [Inf  Inf Inf Inf],...
    'Lower', [0 0 0 -Inf],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
    'StartPoint', [1 round(length(ch1Data)/2) 1 0]...
    );
ft = fittype(model, 'options', fo);
%     xData = (0:length(toFit) - 1)';
[fitobject_ch1, ~, ~] = fit((0:(length(ch1Data) - 1))', ch1Data', ft, fo);
[fitobject_ch2, ~, ~] = fit((0:(length(ch1Data) - 1))', ch2Data', ft, fo);

axes; hold on;
hl = [];
hl(end + 1) = plot(xData, ch1Data, '.', 'Color', mycolors('chat'), 'MarkerSize', 4);
hl(end + 1) = plot(xData, ch2Data, '.', 'Color', mycolors('dat'), 'MarkerSize', 4);
plot(xData, fitobject_ch1.a .* (exp(-1 * ((xData - xData(1)) ./ fitobject_ch1.b) .^fitobject_ch1.c)) + fitobject_ch1.d, '--', 'Color', mycolors('chat'));
plot(xData, fitobject_ch2.a .* (exp(-1 * ((xData - xData(1)) ./ fitobject_ch2.b) .^fitobject_ch2.c)) + fitobject_ch2.d, '--', 'Color', mycolors('dat'));
set(gca, 'YLim', [-4 7]);    
yyaxis right; set(gca, 'YColor', [0 0 0]); ylabel('licks/s');
hl(end + 1) = plot(xData, smoothdata(lickData, 2, 'movmedian', 5),'k');
set(gca, 'YLim', [0 6], 'XLim', xData([1 end]));
yyaxis left;
xlabel('New Cs- trials from reversal')
ylabel('\fontsize{8}Fluor. (\fontsize{12}\sigma\fontsize{8}-baseline)');
% legend({'ACh.', 'Dop.', 'Licks'}, 'Box', 'off', 'Location', 'best', 'Orientation', 'vertical');
formatFigurePublish('size', [1.5 1.2]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end




%% plot the averages, maybe worth showing in supplemental?
smoothWindow = 3;
channels = [1 2];
peakFieldsCh1 = {'phPeakMean_cs_ch1','phPeakMean_us_ch1'};
peakFieldsCh2 = {'phPeakMean_cs_ch2','phPeakMean_us_ch2'};

saveName = [animal '_' strtok(strtok(peakFieldCh1, '_'), '_') '_phCue_revAvg'];
h=ensureFigure(saveName, 1);

trialRange = [-20 40];
% CS MINUS -> CS PLUS
    % x data and indices for baseline (cs- baseline) and ceiling (cs+
    % baseline)
for fieldCounter = 1:2
    peakFieldCh1 = peakFieldsCh1{fieldCounter};
    peakFieldCh2 = peakFieldsCh2{fieldCounter};
    % data fields from RE to plot
    dataFields = {...
        'licks_cs', [0 0 0];...        
        peakFieldCh1, mycolors('chat');...
        peakFieldCh2, mycolors('dat');...
    };
    
    xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
    bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);
    ceilIx(1) = nearest(RE.csMinus.trialsBefore, -30); ceilIx(2) = nearest(RE.csPlus.trialsBefore, 0);
    
    % scale data and plot
    nFields = size(dataFields, 1);
    for counter = 1:nFields
        dataField = dataFields{counter, 1};
        color = dataFields{counter, 2};
        revNorm = [RE.csMinus.(dataField).before RE.csPlus.(dataField).after];
        revNorm = smoothdata(revNorm, 2, 'movmean', smoothWindow, 'omitnan');
    %         baseline = nanmean(RE.csMinus.(dataField).before(:,bl(1):bl(2)), 2); % subtract off Cs- prior to reversal
    %     ceiling = percentile(RE.csPlus.(dataField).before(:,ceilIx(1):ceilIx(2)), 0.8, 2); % divide by Cs+ prior to reversal
    %         revNorm = revNorm - baseline; % RIGHT NOW JUST NORMALIZE, NO BASELINE SUBTRACTION    
    %     revNorm = revNorm ./ ceiling;
        subplot(2, 2, fieldCounter); hold on;
        boundedline(xData, nanmean(revNorm), nanSEM(revNorm), 'cmap', color);
        set(gca, 'XLim', trialRange);        
    end
    %     set(gca, 'XLim', [-50 50]);
    %     xlabel('Trials of new CS+ odor from reversal'); 
    %     ylabel('Cue response normalized'); title([strtok(peakFieldCh1, '_') ', avg ' num2str(nReversals) ' reversals'], 'Interpreter', 'none');

    % CS PLUS -> CS MINUS
        % x data and indices for baseline (cs- baseline) and ceiling (cs+
        % baseline)        
    xData = [RE.csPlus.trialsBefore RE.csMinus.trialsAfter];
    bl(1) = nearest(RE.csMinus.trialsBefore, -30); bl(2) = nearest(RE.csMinus.trialsBefore, 0);
    ceilIx(1) = nearest(RE.csPlus.trialsBefore, -30); ceilIx(2) = nearest(RE.csPlus.trialsBefore, 0);

    % scale data and plot
    for counter = 1:size(dataFields, 1)
        dataField = dataFields{counter, 1};
        color = dataFields{counter, 2};
        revNorm = [RE.csPlus.(dataField).before RE.csMinus.(dataField).after];
        revNorm = smoothdata(revNorm, 2, 'movmean', smoothWindow, 'omitnan');
    %         baseline = nanmean(RE.csMinus.(dataField).before(:,bl(1):bl(2)), 2); % subtract off Cs- prior to reversal
%         ceiling = percentile(RE.csPlus.(dataField).before(:,ceilIx(1):ceilIx(2)), 0.8, 2); % divide by Cs+ prior to reversal
    %         revNorm = revNorm - baseline; % RIGHT NOW JUST NORMALIZE, NO BASELINE SUBTRACTION         
%         revNorm = revNorm ./ ceiling;
        subplot(2,2,2 + fieldCounter); hold on;
        boundedline(xData, nanmean(revNorm), nanSEM(revNorm), 'cmap', color);
        set(gca, 'XLim', trialRange);
    end
end
%     set(gca, 'XLim', [-50 50]);
%     xlabel('Trials of new CS- odor from reversal'); 
%     ylabel('Cue response normalized'); title([strtok(peakFieldCh1, '_') ', avg ' num2str(nReversals) ' reversals'], 'Interpreter', 'none');
    
% 
%     if saveOn
%         saveas(gcf, fullfile(savepath, [saveName '.fig']));
%         saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
%     end


