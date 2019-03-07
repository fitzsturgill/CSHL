saveName = 'cued_reinforcement_allBehavior_v1';
fig = ensureFigure(saveName, 1);
mcLandscapeFigSetup(fig);

%% Define the positions of different axes matrices on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    params = struct();
    matpos_title_1 = [0 0 3/8 .1];
    matpos_title_2 = [3/8 0 (1 - 3/8) .1];
    matpos_rasters = [0 .1 1 0.9];    
    params.cellmargin = [.05 .05 0.05 0.05];    
    params.figmargin = [0.05 0.05 0.025 0.025];
    
%% add titles
params.matpos = matpos_title_1;
[ax txt] = textAxes(fig, 'Cued Reward ------------------------------------', 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
set(txt, 'HorizontalAlignment', 'left', 'Position', [0 .5]);
params.matpos = matpos_title_2;
[ax txt] = textAxes(fig, 'Cued Shock -----------------------------------------------------------------', 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
set(txt, 'HorizontalAlignment', 'left', 'Position', [0 .5]);    

%% cued Reward, cued Shock, all behavior, global numbering
PhotometryField = 'Photometry';
trialNumbering = 'global';
pupField = 'pupDiameterNorm';
CLimFactor = 3;
window = [-4 6];
totalTrials = length(TE.filename);


params.matpos = matpos_rasters;
nRows = 1; nColumns = 9; nAxes = 9;
hax = axesmatrix(nRows, nColumns, 1:nAxes, params, fig);

axes(hax(1)); hold on;
eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
    'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('Licks');         
addStimulusPatch(gca, [0 2]);

% add scatter plot with tags for block type blue = 
c = zeros(totalTrials, 3);
c(TE.BlockNumber == 2, :) = repmat([0 0 1], sum(TE.BlockNumber == 2), 1);
c(TE.BlockNumber == 6, :) = repmat([0.5 0.5 0.5], sum(TE.BlockNumber == 6), 1);
c(TE.BlockNumber == 4, :) = repmat([1 0 0], sum(TE.BlockNumber == 4), 1);
scatter(repmat(window(1) + 0.5, totalTrials, 1), (1:totalTrials)', 4, c);

axes(hax(2)); hold on;
phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('Ch 1'); set(gca, 'XLim', window, 'YTick', []);

if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(3)); hold on;
    phRasterFromTE(TE, trialsByType{1}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('Ch 2'); set(gca, 'XLim', window, 'YTick', []);
else
    axes(hax(3)); hold on;
    alignedDataRaster(TE.pupil.frameAvgNorm, trialsByType{1}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
    title('Eye closure'); set(gca, 'XLim', window, 'YTick', []);
end

axes(hax(4)); hold on;
alignedDataRaster(TE.pupil.(pupField), trialsByType{1}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Pupil diameter'); set(gca, 'XLim', window, 'YTick', []);

axes(hax(5)); hold on;
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('HDB/ChAT'); set(gca, 'XLim', window);

if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(6)); hold on;
    phRasterFromTE(TE, trialsByType{3}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('VTA/DAT'); set(gca, 'XLim', window, 'YTick', []);
end

axes(hax(7)); hold on;
alignedDataRaster(TE.pupil.(pupField), trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Pupil diameter'); set(gca, 'XLim', window, 'YTick', []);
% add scatter plot for shock intensity
try
    c = abs(TE.ShockCurrent);
    scatter(repmat(window(1) + 0.5, totalTrials, 1), (1:totalTrials)', 4, c, 'filled');
catch
end


axes(hax(8)); hold on;
alignedDataRaster(TE.pupil.frameAvgNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Eye closure'); set(gca, 'XLim', window, 'YTick', []);

axes(hax(9)); hold on;
alignedDataRaster(TE.Wheel.data.V, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.Wheel.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.Wheel.Fs, 'CLim', [0 20]);
title('Velocity'); set(gca, 'XLim', window, 'YTick', []);

for counter = 1:length(hax)
    line(hax(counter), [0 0], get(hax(counter), 'YLim'), 'Color', 'r'); line(hax(counter), [2 2], get(hax(counter), 'YLim'), 'Color', 'r');
end
if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveName = 'cued_reinforcement_allBehavior_v2';
fig = ensureFigure(saveName, 1);
mcLandscapeFigSetup(fig);

%% Define the positions of different axes matrices on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    params = struct();
    matpos_title_1 = [0 0 3/8 .1];
    matpos_title_2 = [3/8 0 (1 - 3/8) .1];
    matpos_rasters = [0 .1 1 0.9];    
    params.cellmargin = [.05 .05 0.05 0.05];    
    params.figmargin = [0.05 0.05 0.025 0.025];
    
%% add titles
params.matpos = matpos_title_1;
[ax txt] = textAxes(fig, 'Cued Reward ------------------------------------', 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
set(txt, 'HorizontalAlignment', 'left', 'Position', [0 .5]);
params.matpos = matpos_title_2;
[ax txt] = textAxes(fig, 'Cued Shock -----------------------------------------------------------------', 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
set(txt, 'HorizontalAlignment', 'left', 'Position', [0 .5]);    

%% cued Reward, cued Shock, all behavior, consecutive
PhotometryField = 'Photometry';
trialNumbering = 'consecutive';
CLimFactor = 3;
window = [-4 6];
totalTrials = length(TE.filename);


params.matpos = matpos_rasters;
nRows = 1; nColumns = 8; nAxes = 8;
hax = axesmatrix(nRows, nColumns, 1:nAxes, params, fig);

axes(hax(1)); hold on;
eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
    'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('Licks');         
addStimulusPatch(gca, [0 2]);



axes(hax(2)); hold on;
phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('Ch 1'); set(gca, 'XLim', window, 'YTick', []);



if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(3)); hold on;
    phRasterFromTE(TE, trialsByType{1}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('Ch 2'); set(gca, 'XLim', window, 'YTick', []);
end


axes(hax(4)); hold on;
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('Ch 1'); set(gca, 'XLim', window);

if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(5)); hold on;
    phRasterFromTE(TE, trialsByType{3}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('Ch 2'); set(gca, 'XLim', window, 'YTick', []);
end

axes(hax(6)); hold on;
alignedDataRaster(TE.pupil.pupDiameterNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Pupil diameter'); set(gca, 'XLim', window, 'YTick', []);

axes(hax(7)); hold on;
alignedDataRaster(TE.pupil.frameAvgNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Eye closure'); set(gca, 'XLim', window, 'YTick', []);

axes(hax(8)); hold on;
alignedDataRaster(TE.Wheel.data.V, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.Wheel.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.Wheel.Fs, 'CLim', [0 20]);
title('Velocity'); set(gca, 'XLim', window, 'YTick', []);

for counter = 1:length(hax)
    line(hax(counter), [0 0], get(hax(counter), 'YLim'), 'Color', 'r'); line(hax(counter), [2 2], get(hax(counter), 'YLim'), 'Color', 'r');
end
if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end



%% rasters and averages, combined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveName = 'cued_reinforcement_allBehavior_plusAvgs';
fig = ensureFigure(saveName, 1);
mcLandscapeFigSetup(fig);

%% Define the positions of different axes matrices on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    params = struct();
    matpos_title_1 = [0 0 3/8 .1];
    matpos_title_2 = [3/8 0 (1 - 3/8) .1];
    matpos_rasters = [0 .1 1 0.5];    
    matpos_avgs = [0 0.6 1 0.4];
    params.cellmargin = [.02 .02 0.05 0.05];
    params.figmargin = [0.05 0.05 0.025 0.025];
    
%% add titles
params.matpos = matpos_title_1;
[ax txt] = textAxes(fig, 'Cued Reward -->', 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
set(txt, 'HorizontalAlignment', 'left', 'Position', [0 .5]);
params.matpos = matpos_title_2;
[ax txt] = textAxes(fig, sprintf('       Cued Shock -->    Mouse = %s', subjectName), 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
set(txt, 'HorizontalAlignment', 'left', 'Position', [0 .5]);    

%% cued Reward, cued Shock, all behavior, consecutive
PhotometryField = 'PhotometryExpFit';
trialNumbering = 'consecutive';
CLimFactor = 3;
window = [-4 6];
totalTrials = length(TE.filename);
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);


params.matpos = matpos_rasters;
nRows = 1; nColumns = 9; nAxes = 9;
hax = axesmatrix(nRows, nColumns, 1:nAxes, params, fig);

axes(hax(1)); hold on;
eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
    'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('Licks');         
addStimulusPatch(gca, [0 totalDelay]);



axes(hax(2)); hold on;
phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('Ch 1'); set(gca, 'XLim', window, 'YTick', []);

% add scatter plot with tags for block type blue = 
c = zeros(totalTrials, 3);
c(TE.BlockNumber == 2, :) = repmat([0 0 1], sum(TE.BlockNumber == 2), 1);
c(TE.BlockNumber == 6, :) = repmat([1 1 1], sum(TE.BlockNumber == 6), 1);
c(TE.BlockNumber == 4, :) = repmat([1 0 0], sum(TE.BlockNumber == 4), 1);
% pull out just the rows for trial type 1 (cued Reward) and trial type 2
% (cued Punish)
c_reward = c(trialsByType{1}, :);
c_punish = c(trialsByType{3}, :);
scatter(repmat(window(1) + 0.5, sum(trialsByType{1}), 1), (1:sum(trialsByType{1}))', 4, c_reward, 'filled');

if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(3)); hold on;
    phRasterFromTE(TE, trialsByType{1}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('Ch 2'); set(gca, 'XLim', window, 'YTick', []);
end

axes(hax(4)); hold on;
alignedDataRaster(TE.pupil.(pupField), trialsByType{1}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Pupil diameter'); set(gca, 'XLim', window, 'YTick', []);


axes(hax(5)); hold on;
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
scatter(repmat(window(1) + 0.5, sum(trialsByType{3}), 1), (1:sum(trialsByType{3}))', 4, c_punish, 'filled');
title('Ch 1'); set(gca, 'XLim', window);

if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(6)); hold on;
    phRasterFromTE(TE, trialsByType{3}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('Ch 2'); set(gca, 'XLim', window, 'YTick', []);
end

axes(hax(7)); hold on;
alignedDataRaster(TE.pupil.pupDiameterNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Pupil diameter'); set(gca, 'XLim', window, 'YTick', []);

axes(hax(8)); hold on;
alignedDataRaster(TE.pupil.frameAvgNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Eye closure'); set(gca, 'XLim', window, 'YTick', []);

axes(hax(9)); hold on;
alignedDataRaster(TE.Wheel.data.V, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.Wheel.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.Wheel.Fs, 'CLim', [0 20]);
title('Velocity'); set(gca, 'XLim', window, 'YTick', []);


for counter = 1:length(hax)
    line(hax(counter), [0 0], get(hax(counter), 'YLim'), 'Color', 'r'); line(hax(counter), [totalDelay - 0.1 totalDelay - 0.1], get(hax(counter), 'YLim'), 'Color', 'r');
end



%% averages
% -3 6
params.cellmargin = [.05 .05 0.05 0.05];
window = [-3.5 6];
params.matpos = matpos_avgs;
nRows = 1; nColumns = 3; nAxes = 3;
hax = axesmatrix(nRows, nColumns, 1:nAxes, params, fig);
% for channel = channels

% ch1, aversive
axes(hax(1)); hold on;
[ha, hl] = phPlotAverageFromTE(TE, {trialsByType{[3 6]} trialsByType{4} & TE.BlockNumber == 4}, 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k', 'm'}); %high value, reward
set(gca, 'XLim', window);  ylabel('Ch 1 Fluor. (ZS)'); xlabel('time from cue (s)');
addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
% ch2, aversive
axes(hax(2)); hold on;
try
    [ha, hl] = phPlotAverageFromTE(TE, {trialsByType{[3 6]} trialsByType{4} & TE.BlockNumber == 4}, 2, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k', 'm'}); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    set(gca, 'XLim', window);  ylabel('Ch 2 Fluor. (ZS)');
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
catch
end

% pupil, aversive
xData = TE.pupil.xData + totalDelay; % kludge, because pupil xData is zeroed at Us
pupData = TE.pupil.pupDiameterNorm;
axes(hax(3)); hold on;

boundedline(xData(:), [nanmean(pupData(trialsByType{3}, :)); nanmean(pupData(trialsByType{4} & TE.BlockNumber == 4, :)); nanmean(pupData(trialsByType{6}, :))]',...
    permute([nanstd(pupData(trialsByType{3}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{3},:)), 1));...
    nanstd(pupData(trialsByType{4} & TE.BlockNumber == 4, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{4} & TE.BlockNumber == 4,:)), 1));...
    nanstd(pupData(trialsByType{6}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{6},:)), 1))], [2 3 1]),...
    'cmap', [1 0 0; 0 0 0; 1 0 1]);


 ylabel('Pupil diam.');
set(gca, 'XLim', window);
addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
lh = legend('cued punish', 'omit', 'uncued punish', 'Location', 'best'); set(lh, 'Box', 'off');




if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end





    
