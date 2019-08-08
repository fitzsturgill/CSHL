% FrankenLNL_RewardPunish_exampleMouse

DB = dbLoadExperiment('FrankenLNL_RewardPunish');
animal = 'ACh_3';
savepath = fullfile(DB.path, ['figure' filesep animal]);
ensureDirectory(savepath);
smoothWindow = 1;
saveOn = 1;
window = [-4 4];

photometryField = 'Photometry';
fdField = 'ZS';

success = dbLoadAnimal(DB, animal); % load TE and trial lookups

%%  take a nice swath of sessions for an example figure
sessionNames = {'ACh_3_FrankenLNL_4odors_Jan27_2019_Session1.mat', 'ACh_3_FrankenLNL_4odors_Jan29_2019_Session1.mat', 'ACh_3_FrankenLNL_4odors_Jan30_2019_Session1.mat'};
matches = ismember(TE.filename, sessionNames);
sessionList = unique(TE.sessionIndex(matches));


%% make array of best example traces to choose 1 to show

% rewarding subset

% showTheseR = find(Odor2Valve1Trials & rewardTrials & ismember(TE.sessionIndex, sessionIndex));
% [~, rixR] = sort(rand(size(showTheseR)));
%{
Good rixRs for example traces
sessionIndex = 1
 55    27    49    33    25    14    31     2     5    13
62    51     6    29    15     4    43    31    37    55
65    26     2    49    16    38    14    55    64    35    40    12    41    57    23     7    59    18    37     

sessionIndex = 2
4     1     7    16    26    21    23    34    12     2
 5    21    24    26    22     6     4    16    10    28
 28    12    23    35    26    32     4    22    21    14
  6     2    32    21    25    11    26    15     7    10

sessionIndex = 3
 21     6    12    17     4    31    14    16    18    34
 14    15    12    30     9    34    18    10    21     3
25     5    26     7     4    31    27    18    14     1
%}

exampleArray = [ 55    27    49    33    25    14    31     2     5    13;...
    62    51     6    29    15     4    43    31    37    55;...
    65    26     2    49    16    38    14    55    64    35;...
    4     1     7    16    26    21    23    34    12     2;...
    5    21    24    26    22     6     4    16    10    28;...
    28    12    23    35    26    32     4    22    21    14;...
    6     2    32    21    25    11    26    15     7    10;...
    21     6    12    17     4    31    14    16    18    34;...
    14    15    12    30     9    34    18    10    21     3;...
    25     5    26     7     4    31    27    18    14     1];

exampleArraySessionKey = [1 1 1 2 2 2 2 3 3 3];
    
    


saveName = 'PE_BLA_exampleMouse_traces_array';  
fig = ensureFigure(saveName, 1);    

for counter = 1:length(exampleArraySessionKey)
    showTheseR = find(Odor2Valve1Trials & rewardTrials & ismember(TE.sessionIndex, exampleArraySessionKey(counter)));
    ax = subplot(2,5,counter);
    plot(TE.Photometry.xData - 2, TE.Photometry.data(1).dFF(showTheseR(exampleArray(counter, :)), :)', 'k', 'LineWidth', 0.15); set(gca, 'XLim', window);
    addStimulusPatch(gca, [-2 -1]); addStimulusPatch(gca, [-0.1 0.1]);
    % set(gca, 'YTickLabel', {}); set(gca, 'XTickLabel', {}); 
    set(gca, 'Visible', 'off');
end

formatFigurePublish('size', [1.6 1.1] * 3);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% Choose a good example trace set and plot them for the figure

figSize = [1.6 0.6];
sessionIndex = 2;
whichOne = 7;
rixR = exampleArray(whichOne, :);
showTheseR = find(Odor2Valve1Trials & rewardTrials & ismember(TE.sessionIndex, sessionIndex));

saveName = ['PE_BLA_exampleMouse_traces' num2str(whichOne) '_cued'];  
fig = ensureFigure(saveName, 1);    

ax = subplot(1,1,1);
plot(TE.Photometry.xData - 2, TE.Photometry.data(1).dFF(showTheseR(rixR(1:10)), :)', 'k', 'LineWidth', 0.15); set(gca, 'XLim', window); % xscale shifted relative to reinforcement
addStimulusPatch(gca, [-2 -1]); addStimulusPatch(gca, [-0.1 0.1]);
% set(gca, 'YTickLabel', {}); set(gca, 'XTickLabel', {}); 
set(gca, 'YColor', 'none', 'YTick', [], 'XTickLabel', []);
set(gca, 'YLim', [-0.1 0.3]);

formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end


% also an uncued Example
sessionIndex = 2;
showTheseU = find(uncuedReward & ismember(TE.sessionIndex, sessionIndex));
rixU = [ 7    15    19     1    10     4    18    11     6     8];
saveName = ['PE_BLA_exampleMouse_traces' num2str(whichOne) '_uncued'];  
fig = ensureFigure(saveName, 1);    

ax = subplot(1,1,1);
plot(TE.Photometry.xData - 2, TE.Photometry.data(1).dFF(showTheseU(rixU(1:10)), :)', 'k', 'LineWidth', 0.15); set(gca, 'XLim', window); % xscale shifted relative to reinforcement
addStimulusPatch(gca, [-0.1 0.1]);
% set(gca, 'YTickLabel', {}); set(gca, 'XTickLabel', {}); 
set(gca, 'YColor', 'none', 'YTick', [], 'XTickLabel', []);
set(gca, 'YLim', [-0.1 0.3]);

formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% cued reward rasters, cued and uncued for photometry



fdField = 'ZS';
tcolor = mycolors('chat');
sessionIndices = [2 3];
hitTrials = TE.licks_cs.rate > 0;
figSize = [1.6 0.81];
saveName = 'PE_BLA_exampleMouse_Appetitive_Rasters_cuedLicks';  
ensureFigure(saveName, 1);
subplot(1,1,1);
[~, lh] = eventRasterFromTE(TE, trialsByType{1} & hitTrials & ismember(TE.sessionIndex, sessionIndices), 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
set(gca, 'XLim', window);
set(gca, 'YTick', [1 50]);
% ylabel('Trial #');
formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end        


saveName = 'PE_BLA_exampleMouse_Appetitive_Rasters_cued';  
ensureFigure(saveName, 1);
climfactor = 3;  
subplot(1,1,1); phRasterFromTE(TE, trialsByType{1} & hitTrials & ismember(TE.sessionIndex, sessionIndices), 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry', 'zeroTimes', TE.Us, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
set(gca, 'YTick', [1 50]);
formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end        

saveName = 'PE_BLA_exampleMouse_Appetitive_Rasters_uncued';  
ensureFigure(saveName, 1);
subplot(1,1,1); phRasterFromTE(TE, uncuedReward & hitTrials & ismember(TE.sessionIndex, sessionIndices), 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'Photometry', 'zeroTimes', TE.Us, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
set(gca, 'YTick', [1 30]);
formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end        
%% averages, appetitive, one side only (left side, ch = 1)
figSize = [1.7 0.9];
fdField = 'ZS';
saveName = 'PE_BLA_exampleMouse_avgs_simple';  
h=ensureFigure(saveName, 1); 
% use the same sessions used for the rasters
sessionIndices = [2 3];
linecolors = [0 0 1; 0 0 0; 0 1 1];            

subplot(1, 1, 1);
% set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');

[ha, hl] = phPlotAverageFromTE(TE, {trialsByType{1} & ismember(TE.sessionIndex, sessionIndices), trialsByType{2} & ismember(TE.sessionIndex, sessionIndices), trialsByType{5} & ismember(TE.sessionIndex, sessionIndices)}, 1,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window, 'cmap', linecolors, 'alpha', 1); %high value, reward
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'best'); legend('boxoff');
ylabel('F(\fontsize{10}\sigma\fontsize{7}-baseline)');  set(gca, 'XLim', window);
xlabel('Time from reinforcement (s)');

formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
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
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
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
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
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
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end


return
%% Development code below:
%% normalize responses to punishment by those to reward and compare between left and right BLA

saveName = 'LeftRight_PunNorm_dev';
ensureFigure(saveName, 1);
formatFigurePublish('size', [4 2], 'fontSize', 8);
animals = {'ACh_7', 'ACh_3'};


% xlims = [min(TE.phPeakMean_cs(2).data(allTrials)) max(TE.phPeakMean_cs(2).data(allTrials)); min(TE.phPeakMean_us(2).data(allTrials)) max(TE.phPeakMean_us(2).data(allTrials))];

for acounter = 1:length(animals)
    subplot(1,2,acounter); hold on;
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
        xData = TE.phPeakMean_us(1).data(trialSets(:, counter)); yData = TE.phPeakMean_us(2).data(trialSets(:, counter)); 
        xData = xData(:); yData = yData(:);
        if counter == 1
            xDenom = nanmean(xData);
            yDenom = nanmean(yData);
        end
        xData = xData ./ xDenom;
        yData = yData ./ xDenom;

        h(end + 1) = scatter(xData, yData, 20, linecolors(counter, :), '.', 'MarkerFaceColor', 'flat', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);

        xDataStruct(counter).Mean  = nanmean(xData);
        xDataStruct(counter).SEM = nanstd(xData) / sqrt(sum(isfinite(xData)));
        yDataStruct(counter).Mean = nanmean(yData);
        yDataStruct(counter).SEM = nanstd(yData) / sqrt(sum(isfinite(yData)));

    end
    h = [];
    for counter = 1:size(trialSets, 2)
        errorbar(xDataStruct(counter).Mean, yDataStruct(counter).Mean, -yDataStruct(counter).SEM, yDataStruct(counter).SEM,-xDataStruct(counter).SEM, xDataStruct(counter).SEM, 'Color', linecolors(counter, :), 'LineWidth', 1, 'CapSize', 2, 'Marker', 'none');
        h(end + 1) = plot(xDataStruct(counter).Mean, yDataStruct(counter).Mean, '-', 'Color', linecolors(counter, :));
    end    
    
        sameXYScale(gca);
        addOrginLines(gca);
        legend(h, trialSetNames, 'Location', 'Best'); 
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
    export_fig(fullfile(savepath, saveName), '-eps');
end

%% dev, fano factor, skew, kurtosis
sessionIndices = [2 3];
hitTrials = TE.licks_cs.rate > 0;

[cuedData, xData] = phAlignedWindow(TE, trialsByType{1} & hitTrials & ismember(TE.sessionIndex, sessionIndices), 1,...
    'zeroField', 'Us', 'FluorDataField', 'raw', 'PhotometryField', 'Photometry', 'zeroTimes', TE.Us, 'window', [-2 2]);


data_std = std(cuedData, 0, 1);
data_mean = mean(cuedData);
data_skew = skewness(cuedData, 1);

data_cv = data_std ./ data_mean;
data_fano = data_std.^2 ./ data_mean;

ensureFigure('test', 1); 
subplot(2,2,1); plot(xData, data_mean);
subplot(2,2,2); plot(xData, data_std);
subplot(2,2,3); plot(xData, data_skew);
subplot(2,2,4); plot(xData, data_cv);

%%
ensureFigure('test2', 1);
axes;
yyaxis left;
plot(xData, data_mean);
yyaxis right;
plot(xData, data_cv);
grid on;

%%
[cuedData, xData] = phAlignedWindow(TE, uncuedReward & hitTrials & ismember(TE.sessionIndex, sessionIndices), 1,...
    'zeroField', 'Us', 'FluorDataField', 'raw', 'PhotometryField', 'Photometry', 'zeroTimes', TE.Us, 'window', window);


data_std = std(cuedData, 0, 1);
data_mean = mean(cuedData);
data_skew = skewness(cuedData, 1);
data_cv = data_std ./ data_mean;
data_fano = data_std.^2;

ensureFigure('test', 1); 
subplot(2,2,1); plot(xData, data_mean);
subplot(2,2,2); plot(xData, data_std);
subplot(2,2,3); plot(xData, data_skew);
subplot(2,2,4); plot(xData, data_cv);

%%
ensureFigure('test2', 1);
axes;
yyaxis left;
plot(xData, data_mean);
yyaxis right;
plot(xData, data_fano);
grid on;




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
