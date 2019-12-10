% FrankenLNL_RewardPunish_exampleMouse

DB = dbLoadExperiment('FrankenLNL_RewardPunish');
animal = 'ACh_3';
savepath = fullfile(DB.path, ['figure' filesep animal]);
ensureDirectory(savepath);
smoothWindow = 1;
saveOn = 1;
window = [-3 3];

photometryField = 'PhotometryHF';
fdField = 'ZS';

success = dbLoadAnimal(DB, animal); % load TE and trial lookups

%%  take a nice swath of sessions for an example figure
sessionNames = {'ACh_3_FrankenLNL_4odors_Jan27_2019_Session1.mat', 'ACh_3_FrankenLNL_4odors_Jan29_2019_Session1.mat', 'ACh_3_FrankenLNL_4odors_Jan30_2019_Session1.mat'};
matches = ismember(TE.filename, sessionNames);
sessionList = unique(TE.sessionIndex(matches));




%% cued reward rasters, cued and uncued for photometry

fdField = 'ZS';
tcolor = mycolors('chat');
sessionIndices = [2 3];
hitTrials = TE.licks_cs.rate > 0;
figSize = [1.6 0.81];
nTrials = 30;


saveName = 'Fig1_Sensor_cuedRasters';  
ensureFigure(saveName, 1);
climfactor = 3;  
trials = find(trialsByType{1} & hitTrials & ismember(TE.sessionIndex, sessionIndices));
subplot(1,1,1); phRasterFromTE(TE, trials(1:nTrials), 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'zeroTimes', TE.Us, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
set(gca, 'YTick', [nTrials], 'XLim', window, 'XTick', [-3 0 3], 'XTickLabel', {});
formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    export_fig(fullfile(savepath, saveName), '-eps');
end        

saveName = 'Fig1_Sensor_uncuedRasters';  
ensureFigure(saveName, 1);
trials = find(uncuedReward & hitTrials & ismember(TE.sessionIndex, sessionIndices));
subplot(1,1,1); phRasterFromTE(TE, trials(1:nTrials), 1, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'zeroTimes', TE.Us, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
set(gca, 'YTick', [nTrials], 'XLim', window, 'XTick', [-3 0 3], 'XTickLabel', {});
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
set(gca, 'XLim', window, 'XTick', [-3 0 3]);
addStimulusPatch(gca, [-2 -1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'best'); legend('boxoff');
% ylabel('F(\fontsize{10}\sigma\fontsize{7}-baseline)');  set(gca, 'XLim', window);
% xlabel('Time from reinforcement (s)');

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
