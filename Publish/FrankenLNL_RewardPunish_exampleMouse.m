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
figSize = [0.7 0.4];
nTraces = 1;
window = [-3 7];
showTheseR = find(Odor2Valve1Trials & rewardTrials & ismember(TE.sessionIndex, sessionIndex));
[~, rixR] = sort(rand(size(showTheseR)));

saveName = 'PE_BLA_exampleMouse_traces_LR';  
fig = ensureFigure(saveName, 1);    
ax = subplot(1,1,1); hold on;
% trial 16 is a nice one for the figure
rixR(1) = 16; % for now, use 16
ch1Data = TE.Photometry.data(1).dFF(showTheseR(rixR(1:nTraces)), :)';
ch2Data = TE.Photometry.data(2).dFF(showTheseR(rixR(1:nTraces)), :)' + range(ch1Data);
set(gca, 'YLim', [min(ch1Data) - range(ch1Data) * 0.1 max(ch2Data) + range(ch1Data) * 0.1]);
hp = addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 1); 
hp = addStimulusPatch(gca, [2.9 3.1], '', [0.8 0.8 0.8], 1);
plot(TE.Photometry.xData(21:end), ch1Data(21:end,:), 'k', 'LineWidth', 0.4); set(gca, 'XLim', window);
plot(TE.Photometry.xData(21:end), ch2Data(21:end,:), 'k', 'LineWidth', 0.4); set(gca, 'XLim', window);

set(gca, 'Visible', 'off');

formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg'])); 
end



%% 
saveName = 'PE_BLA_exampleMouse_traces_uncued';  
fig = ensureFigure(saveName, 1);    

% aversive subset
showTheseU = find(uncuedReward & ismember(TE.sessionIndex, [2]));
[~, rixU] = sort(rand(size(showTheseU)));

ax = subplot(1,1,1);
plot(TE.Photometry.xData, TE.Photometry.data(1).dFF(showTheseU(rixU(1:10)), :)', 'k', 'LineWidth', 0.15); set(gca, 'XLim', [-2 6]);
 addStimulusPatch(gca, [1.9 2.1]);
% set(gca, 'YTickLabel', {}); set(gca, 'XTickLabel', {}); 
set(gca, 'Visible', 'off');
set(gca, 'YLim', [-0.1 0.3]);

% formatFigurePublish('size', [1.6 1.1] );

% if saveOn 
%     export_fig(fullfile(savepath, saveName), '-eps');
% end

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
FluorField_cs = 'phPeakMean_cs';
FluorField_us = 'phPeakMean_us';
% xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];

subplot(1,2,1); hold on; subplot(1,2,2); hold on;
for counter = 1:size(trialSets, 2)
    h=[];
    h(1) = subplot(1,2,1); %set(gca, 'XLim', xlims(1,:));
%     allTrials = allTrials | trialSets{counter};
    xData = TE.(FluorField_cs)(1).data(trialSets(:, counter)); yData = TE.(FluorField_cs)(2).data(trialSets(:, counter)); 
%     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for cs
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off; %,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));

    h(2) = subplot(1,2,2); %set(gca, 'XLim', xlims(2,:));
    xData = TE.(FluorField_us)(1).data(trialSets(:, counter)); yData = TE.(FluorField_us)(2).data(trialSets(:, counter)); 
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
    xData = TE.(FluorField_cs)(1).data(trialSets(:, counter)) - nanmean(TE.(FluorField_cs)(1).data(trialSets(:, counter))); 
    yData = TE.(FluorField_cs)(2).data(trialSets(:, counter)) - nanmean(TE.(FluorField_cs)(2).data(trialSets(:, counter))); 
%     scatter(TE.phPeakMean_cs(2).data(trialSets{counter}), TE.phPeakMean_cs(1).data(trialSets{counter}), 8, linecolors(counter, :), '.');
    scatter(xData, yData, 8, linecolors(counter, :), '.');
    % fit for cs
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData, yData, 'poly1', fo); 
    fph=plot(fob); legend off; %,'predfunc'); legend off;
    set(fph, 'LineWidth', 0.5, 'Color', linecolors(counter, :));

    h(2) = subplot(1,2,2); %set(gca, 'XLim', xlims(2,:));
    xData = TE.(FluorField_us)(1).data(trialSets(:, counter)) - nanmean(TE.(FluorField_us)(1).data(trialSets(:, counter))); 
    yData = TE.(FluorField_us)(2).data(trialSets(:, counter)) - nanmean(TE.(FluorField_us)(2).data(trialSets(:, counter))); 
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
%% Development code below: normalize responses to punishment by those to reward and compare between left and right BLA


FluorField_us = 'phPeakMean_us';

saveName = sprintf('LeftRight_PunNorm_%s', FluorField_us);
ensureFigure(saveName, 1);
% formatFigurePublish('size', [4 2], 'fontSize', 8);
% animals = {'ACh_7', 'ACh_3'};
animals = DB.animals;


% xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];

nCol = ceil(sqrt(length(animals)));
nRow = ceil(length(animals) / nCol);

for acounter = 1:length(animals)
    subplot(nCol, nRow,acounter); hold on;
    success = dbLoadAnimal(DB, animals{acounter}); % load TE and trial lookups
    title(animals{acounter}, 'Interpreter', 'none');
    % reward is first in the trialSets list, use it to normalize
    linecolors = [0 0 1; 1 0 0; mycolors('shock')];         
    trialSets = [uncuedReward, uncuedPunish, uncuedShock];
    trialSetNames = {'Reward', 'Air Puff', 'Shock'};
    allTrials = sum(trialSets, 2) ~= 0;    
    xDataStruct = struct('Mean', cell(size(trialSets, 2), 1), 'SEM', cell(size(trialSets, 2), 1));
    yDataStruct = struct('Mean', cell(size(trialSets, 2), 1), 'SEM', cell(size(trialSets, 2), 1));
    h=[];   
    for counter = 1:size(trialSets, 2)        
    %     allTrials = allTrials | trialSets{counter};
        xData = TE.(FluorField_us)(1).data(trialSets(:, counter)); yData = TE.(FluorField_us)(2).data(trialSets(:, counter)); 
        xData = xData(:); yData = yData(:);
        if counter == 1
            xDenom = nanmean(xData);
            yDenom = nanmean(yData);
        end
        xData = xData ./ xDenom;
        yData = yData ./ yDenom;

        h(end + 1) = scatter(xData, yData, 20, linecolors(counter, :), '.', 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);

        xDataStruct(counter).Mean  = nanmean(xData);
        xDataStruct(counter).SEM = nanstd(xData) / sqrt(sum(isfinite(xData)));
        yDataStruct(counter).Mean = nanmean(yData);
        yDataStruct(counter).SEM = nanstd(yData) / sqrt(sum(isfinite(yData)));
    end
    addUnityLine(gca, [0.3 0.3 0.3]);

%     xlim = get(gca,'XLim');
%     ylim = get(gca,'YLim');
% 
%     p1 = min([min(xlim) min(ylim)]);
%     p2 = min([max(xlim) max(ylim)]);
%     p1 = 

%     h = line('Parent',gca,'XData',[p1 p2],'YData',[p2 p1]);

%     set(h,'Color',[0.3 0.3 0.3]);
    h = [];
    for counter = 1:size(trialSets, 2)
        errorbar(xDataStruct(counter).Mean, yDataStruct(counter).Mean, -yDataStruct(counter).SEM, yDataStruct(counter).SEM,-xDataStruct(counter).SEM, xDataStruct(counter).SEM, 'Color', linecolors(counter, :), 'LineWidth', 1, 'CapSize', 2, 'Marker', 'none');
        h(end + 1) = plot(xDataStruct(counter).Mean, yDataStruct(counter).Mean, '-', 'Color', linecolors(counter, :));
    end
    
    if acounter == 2
        set(gca, 'YLim', [-5 6], 'XLim', [-5 6]);
    end
    sameXYScale(gca);
    addOrginLines(gca);
    legend(h, trialSetNames, 'Location', 'Best'); legend('boxoff'); 
    % subplot(1,2,1);
    % textBox('Cue', [], [0.5 0.95], 8);
    % xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
%     ylabel('\fontsize{8}Right (\fontsize{12}\sigma\fontsize{8}-baseline)');
    ylabel('Right (reward-normalized)');
    % subplot(1,2,2);
    % textBox('Outcome', [], [0.5 0.95], 8);
%     xlabel('\fontsize{8}Left (\fontsize{12}\sigma\fontsize{8}-baseline)');
    xlabel('Left (reward-normalized)');

    % ylabel('');   
end





if saveOn 
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));

    export_fig(gcf,fullfile(savepath, [saveName '.pdf']));  % write to pdf
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
