% FrankenLNL_RewardPunish_exampleMouse

DB = dbLoadExperiment('FrankenLNL_RewardPunish');
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);
smoothWindow = 1;
saveOn = 1;
animal = 'ACh_7';

photometryField = 'Photometry';
fdField = 'ZS';

success = dbLoadAnimal(DB, animal); % load TE and trial lookups

%% % find the session index for example session with air puff
sessionName = 'ACh_7_FrankenLNL_4odors_Feb27_2019_Session1.mat';
matches = strcmp(TE.filename, sessionName);
sessionIndex = unique(TE.sessionIndex(matches));


%% example traces
% rewarding subset
window = [-4 7];
showTheseR = find(Odor2Valve1Trials & rewardTrials & TE.sessionIndex == sessionIndex);
[~, rixR] = sort(rand(size(showTheseR)));



saveName = 'PE_BLA_exampleMouse_traces';  
fig = ensureFigure(saveName, 1);    


ax = subplot(1,1,1);
plot(TE.Photometry.xData, TE.Photometry.data(1).dFF(showTheseR(rixR(1:10)), :)', 'k', 'LineWidth', 0.15); set(gca, 'XLim', window);
addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [2.9 3.1]);
% set(gca, 'YTickLabel', {}); set(gca, 'XTickLabel', {}); 
set(gca, 'Visible', 'off');
% title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Left']);


% 
% ax(2) = subplot(1,2,2);
% plot(TE.Photometry.xData, TE.Photometry.data(2).dFF(showTheseR(rixR(1:10)), :)', 'k', 'LineWidth', 0.15); set(gca, 'XLim', window);
% addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [2.9 3.1]);
% % title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Left']);
% % set(gca, 'YTickLabel', {}); set(gca, 'XTickLabel', {}); 
% set(gca, 'Visible', 'off');


% aversive subset
showTheseA = find(Odor2Valve2Trials & punishTrials & TE.sessionIndex == sessionIndex);
[~, rixA] = sort(rand(size(showTheseA)));

formatFigurePublish('size', [1.6 1.1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% cued reward rasters, sorted by first lick
% computer latency to first lick


fdField = 'ZS';
tcolor = mycolors('chat');
window = [-4 7]; % relative to Us
TE.firstLick = calcEventLatency(TE, 'Port1In', TE.Cue2, TE.Us); % my calcEventLatency functions computes the latency between Bpod time stamps and events (e.g. licks and the start of a bpod state)
hitTrials = TE.licks_cs.rate > 0;
saveName = 'PE_BLA_exampleMouse_cuedReward_Rasters_sorted';  
ensureFigure(saveName, 1);

subplot(1,3,1);
[~, lh] = eventRasterFromTE(TE, trialsByType{1} & hitTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroField', 'Cue2', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.firstLick);
set(gca, 'XLim', window);
% set(lh, 'LineWidth', 0.3, 'Color', [0 0 0]);
title('licking');
%     xlabel('Time from odor (s)');
climfactor = 3;  

lickOnsets = TE.firstLick(trialsByType{1});
lickOnsets = sort(lickOnsets);
subplot(1,3,2); phRasterFromTE(TE, trialsByType{1} & hitTrials, 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry', 'sortValues', TE.firstLick, 'zeroTimes', TE.Cue2, 'window', window); % 'CLimFactor', CLimFactor,
line(lickOnsets, (1:sum(trialsByType{1}))', 'Parent', gca, 'Color', 'r', 'LineWidth', 2);
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Left']);
subplot(1,3,3); phRasterFromTE(TE, trialsByType{1} & hitTrials, 2, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry', 'sortValues', TE.firstLick, 'zeroTimes', TE.Cue2, 'window', window); % 'CLimFactor', CLimFactor,
line(lickOnsets, (1:sum(trialsByType{1}))', 'Parent', gca, 'Color', 'r', 'LineWidth', 2);
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Right']);


        
%% averages, appetitive

fdField = 'ZS';
saveName = 'PE_BLA_exampleMouse_avgs';  
h=ensureFigure(saveName, 1); 
sessionIndexList = 5; % just the 1 session because early sessions surprise modulation hasn't developed whereas later sessions I shorten the delay (surprise modulation is consistent) but it messes up the graph having 2 different delays
window = [-7 4];
linecolors = [0 0 1; 0 0 0; 0 1 1];            

subplot(1, 2, 1);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Left']);
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{2} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 1,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window, 'cmap', linecolors); %high value, reward

% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northwest'); legend('boxoff');
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', window);

subplot(1, 2, 2);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Right']);
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{2} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{5} & ismember(TE.sessionIndex, sessionIndexList)}, 2,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window, 'cmap', linecolors); %high value, reward

% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northwest'); legend('boxoff');
set(gca, 'XLim', window);


formatFigurePublish('size', [3 1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% averages, aversive, airpuff, also add uncued shock
fdField = 'ZS';
saveName = 'PE_BLA_exampleMouse_aversive_avgs';  
h=ensureFigure(saveName, 1); 
sessionIndexList = 5;


linecolors = [1 0 0; 0 0 0; 1 0 1];            
window = [-7 4];
subplot(1, 2, 1);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {trialsByType{3} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{4} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{6} & ismember(TE.sessionIndex, sessionIndexList)}, 1,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window, 'cmap', linecolors); hold on; %high value, reward
% add uncued shock
[ha, hl] = phPlotAverageFromTE(TE, uncuedShock, 1,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', [-5 4], 'cmap', mycolors('shock')); 

% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northwest'); legend('boxoff');
ylabel('F(\fontsize{12}\sigma\fontsize{8}-baseline)'); xlabel('Time from');  set(gca, 'XLim', window);

subplot(1, 2, 2);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {trialsByType{3} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{4} & ismember(TE.sessionIndex, sessionIndexList), trialsByType{6} & ismember(TE.sessionIndex, sessionIndexList)}, 2,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window, 'cmap', linecolors); hold on; %high value, reward
% add uncued shock
[ha, hl] = phPlotAverageFromTE(TE, uncuedShock, 2,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', [-5 4], 'cmap', mycolors('shock')); 

% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northwest'); legend('boxoff');
xlabel('reinforcement (s)'); set(gca, 'XLim', window);


formatFigurePublish('size', [3 1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% SIGNAL CORRELATIONS: show that cue and  and reward responses are correlated on a trial-by-trial basis between let and right BLA

saveName = 'LeftRight_BLA_SIGNAL_correlations';
ensureFigure(saveName, 1);
hitTrials = TE.licks_cs.rate > 0;
linecolors = [0 0 1; 0 1 1; 1 0 0; 0 1 0; mycolors('shock')];         
trialSets = [Odor2Valve1Trials & hitTrials & rewardTrials, uncuedReward, Odor2Valve2Trials & punishTrials, Odor2Valve2Trials & shockTrials];
allTrials = sum(trialSets, 2) ~= 0;
% xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];

subplot(1,2,1); hold on; subplot(1,2,2); hold on;
for counter = 1:size(trialSets, 2)
    h=[];
    h(1) = subplot(1,2,1); %set(gca, 'XLim', xlims(1,:));
%     allTrials = allTrials | trialSets{counter};
    xData = TE.phPeakMean_cs(1).data(trialSets(:, counter)); yData = TE.phPeakMean_cs(2).data(trialSets(:, counter)); 
%     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for cs
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off; %,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));

    h(2) = subplot(1,2,2); %set(gca, 'XLim', xlims(2,:));
    xData = TE.phPeakMean_us(1).data(trialSets(:, counter)); yData = TE.phPeakMean_us(2).data(trialSets(:, counter)); 
%     scatter(TE.phPeakMean_us(2).data(trialSets{counter}), TE.phPeakMean_us(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');    
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for us
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off;% ,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
    sameXYScale(h);
end
subplot(1,2,1);
textBox('Cue', [], [0.5 0.95], 8);
xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
ylabel('\fontsize{8}Right (\fontsize{12}\sigma\fontsize{8}-baseline)');
addOrginLines;
subplot(1,2,2);
textBox('Outcome', [], [0.5 0.95], 8);
xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
ylabel('');
addOrginLines;

formatFigurePublish('size', [2.5 1.1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end


%% NOISE CORRELATIONS: show that cue and reward responses are correlated on a trial-by-trial basis between let and right BLA

saveName = 'LeftRight_BLA_NOISE_correlations';
ensureFigure(saveName, 1);
hitTrials = TE.licks_cs.rate > 0;
linecolors = [0 0 1; 0 1 1; 1 0 0; 0 1 0; mycolors('shock')];         
trialSets = [Odor2Valve1Trials & hitTrials & rewardTrials, uncuedReward, Odor2Valve2Trials & punishTrials, Odor2Valve2Trials & shockTrials];
allTrials = sum(trialSets, 2) ~= 0;
% xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];

subplot(1,2,1); hold on; subplot(1,2,2); hold on;
for counter = 1:size(trialSets, 2)
    h=[];
    h(1) = subplot(1,2,1); %set(gca, 'XLim', xlims(1,:));
%     allTrials = allTrials | trialSets{counter};
    xData = TE.phPeakMean_cs(1).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_cs(1).data(trialSets(:, counter))); 
    yData = TE.phPeakMean_cs(2).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_cs(2).data(trialSets(:, counter))); 
%     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for cs
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off; %,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));

    h(2) = subplot(1,2,2); %set(gca, 'XLim', xlims(2,:));
    xData = TE.phPeakMean_us(1).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_us(1).data(trialSets(:, counter))); 
    yData = TE.phPeakMean_us(2).data(trialSets(:, counter)) - nanmean(TE.phPeakMean_us(2).data(trialSets(:, counter))); 
%     scatter(TE.phPeakMean_us(2).data(trialSets{counter}), TE.phPeakMean_us(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');    
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for us
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off;% ,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
    sameXYScale(h);
end
subplot(1,2,1);
textBox('Cue', [], [0.5 0.95], 8);
xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
ylabel('\fontsize{8}Right (\fontsize{12}\sigma\fontsize{8}-baseline)');
subplot(1,2,2);
textBox('Outcome', [], [0.5 0.95], 8);
xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
ylabel('');

formatFigurePublish('size', [2.5 1.1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end


return
%% Development code below:
%% normalize responses to punishment by those to reward and compare between left and right BLA

saveName = 'LeftRight_PunNorm_dev';
ensureFigure(saveName, 1);


% reward is first in the trialSets list, use it to normalize
linecolors = [0 0 1; 1 0 0; mycolors('shock')];         
trialSets = [uncuedReward, uncuedPunish, uncuedShock];
allTrials = sum(trialSets, 2) ~= 0;
% xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];

subplot(1,1,1); hold on; 
for counter = 1:size(trialSets, 2)
    h=[];   
%     allTrials = allTrials | trialSets{counter};
    xData = TE.phPeakMean_us(1).data(trialSets(:, counter)); yData = TE.phPeakMean_us(2).data(trialSets(:, counter)); 
    xData = xData(:); yData = yData(:);
    if counter == 1
        xDenom = nanmean(xData);
        yDenom = nanmean(yData);
    end
    xData = xData ./ xDenom;
    yData = yData ./ xDenom;
    scatter(xData, yData, 20, linecolors(counter, :), '.', 'MarkerFaceColor', 'flat');

    xMean = nanmean(xData);
    xSEM = nanstd(xData) / sqrt(sum(isfinite(xData)));
    yMean = nanmean(yData);
    ySEM = nanstd(yData) / sqrt(sum(isfinite(yData)));

    errorbar(xMean, yMean, -ySEM, ySEM,...
        -xSEM, xSEM, 'Color', linecolors(counter, :), 'LineWidth', 2);
%     errorbar(xMean, yMean, -0.5, 0.5,...
%         -0.5, 0.5, 'o', 'Color', linecolors(counter, :));


    
    % fit for us
%     fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
%     fob = fit(xData, yData, 'poly1', fo); 
%     fph=plot(fob, 'predfunc'); legend off;
%     set(fph, 'LineWidth', 1, 'Color', linecolors(counter, :));

%     h(2) = subplot(1,2,2); %set(gca, 'XLim', xlims(2,:));
%     xData = TE.phPeakMean_us(1).data(trialSets(:, counter)); yData = TE.phPeakMean_us(2).data(trialSets(:, counter)); 
% %     scatter(TE.phPeakMean_us(2).data(trialSets{counter}), TE.phPeakMean_us(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');    
%     scatter(xData, yData, 8, linecolors(counter, :), '.');
%     % fit for us
%     fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
%     fob = fit(xData, yData, 'poly1', fo); 
%     fph=plot(fob); legend off;% ,'predfunc'); legend off;
%     set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));

end
    sameXYScale(gca);
    addOrginLines(gca);
% subplot(1,2,1);
% textBox('Cue', [], [0.5 0.95], 8);
% xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
ylabel('\fontsize{8}Right (\fontsize{12}\sigma\fontsize{8}-baseline)');
% subplot(1,2,2);
% textBox('Outcome', [], [0.5 0.95], 8);
xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
% ylabel('');

formatFigurePublish('size', [2 2]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end



%% scrap:

% %% SIGNAL CORRELATIONS II: ACTUAL signal correlations (doesn't really look good)
% 
% saveName = 'LeftRight_BLA_SIGNAL_correlations_II';
% ensureFigure(saveName, 1);
% hitTrials = TE.licks_cs.rate > 0;
% linecolors = [0 0 1; 0 1 1; 1 0 0; 0 1 0;];%mycolors('shock')];
% trialSets = [Odor2Valve1Trials & hitTrials & rewardTrials, uncuedReward, Odor2Valve2Trials & punishTrials, Odor2Valve2Trials & shockTrials];
% allTrials = sum(trialSets, 2) ~= 0;
% % xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];
% 
% 
% % gather the data
% nConditions = size(trialSets, 2); 
% 
% xCue = struct(...
%     'avg', zeros(nConditions, 1),...
%     'sem', zeros(nConditions, 1)...
%     );
% yCue = xCue;
% xUs = xCue;
% yUs = yCue;
% 
% for counter = 1:size(trialSets, 2)
%     h=[];    
% %     allTrials = allTrials | trialSets{counter};
%     xCue.avg(counter) = nanmean(TE.phPeakMean_cs(1).data(trialSets(:, counter)));
%     yCue.avg(counter) = nanmean(TE.phPeakMean_cs(2).data(trialSets(:, counter)));
% 
%     xUs.avg(counter) = nanmean(TE.phPeakMean_us(1).data(trialSets(:, counter)));
%     yUs.avg(counter) = nanmean(TE.phPeakMean_us(2).data(trialSets(:, counter)));
% 
% % %     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
% %     scatter(xData, yData, 8, linecolors(counter, :), '.');
% %     % fit for cs
% %     fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
% %     fob = fit(xData, yData, 'poly1', fo); 
% %     fph=plot(fob); legend off; %,'predfunc'); legend off;
% %     set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
% % 
% %     h(2) = subplot(1,2,2); %set(gca, 'XLim', xlims(2,:));
% %     xData = TE.phPeakMean_us(1).data(trialSets(:, counter)); yData = TE.phPeakMean_us(2).data(trialSets(:, counter)); 
% % %     scatter(TE.phPeakMean_us(2).data(trialSets{counter}), TE.phPeakMean_us(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');    
% %     scatter(xData, yData, 8, linecolors(counter, :), '.');
% %     % fit for us
% %     fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
% %     fob = fit(xData, yData, 'poly1', fo); 
% %     fph=plot(fob); legend off;% ,'predfunc'); legend off;
% %     set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));
% %     sameXYScale(h);
% end
% 
% subplot(1,2,1); hold on; 
% scatter(xCue.avg, yCue.avg, 30, linecolors, 'o', 'MarkerFaceColor', 'flat');
% subplot(1,2,2); hold on;
% scatter(xUs.avg, yUs.avg, 30, linecolors, 'o', 'MarkerFaceColor', 'flat');
% 
% 
% % subplot(1,2,1);
% % textBox('Cue', [], [0.5 0.95], 8);
% % xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
% % ylabel('\fontsize{8}Right (\fontsize{12}\sigma\fontsize{8}-baseline)');
% % subplot(1,2,2);
% % textBox('Outcome', [], [0.5 0.95], 8);
% % xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
% % ylabel('');
% % 
% % formatFigurePublish('size', [2.5 1.1]);
% % 
% % if saveOn 
% %     export_fig(fullfile(savepath, saveName), '-eps');
% % end
