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

% add scatter plot with tags for block type blue = 
c = zeros(totalTrials, 3);
c(TE.BlockNumber == 2, :) = repmat([0 0 1], sum(TE.BlockNumber == 2), 1);
c(TE.BlockNumber == 6, :) = repmat([0.5 0.5 0.5], sum(TE.BlockNumber == 6), 1);
c(TE.BlockNumber == 4, :) = repmat([1 0 0], sum(TE.BlockNumber == 4), 1);
scatter(repmat(window(1) + 0.5, totalTrials, 1), (1:totalTrials)', 4, c);

axes(hax(2)); hold on;
phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('HDB/ChAT'); set(gca, 'XLim', window, 'YTick', []);

if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(3)); hold on;
    phRasterFromTE(TE, trialsByType{1}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('VTA/DAT'); set(gca, 'XLim', window, 'YTick', []);
else
    axes(hax(3)); hold on;
    alignedDataRaster(TE.pupil.frameAvgNorm, trialsByType{1}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
    title('Eye closure'); set(gca, 'XLim', window, 'YTick', []);
end


axes(hax(4)); hold on;
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('HDB/ChAT'); set(gca, 'XLim', window);

if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(5)); hold on;
    phRasterFromTE(TE, trialsByType{3}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('VTA/DAT'); set(gca, 'XLim', window, 'YTick', []);
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

% add scatter plot with tags for block type blue = 
% c = zeros(totalTrials, 3);
% c(TE.BlockNumber == 2, :) = repmat([0 0 1], sum(TE.BlockNumber == 2), 1);
% c(TE.BlockNumber == 6, :) = repmat([0.5 0.5 0.5], sum(TE.BlockNumber == 6), 1);
% c(TE.BlockNumber == 4, :) = repmat([1 0 0], sum(TE.BlockNumber == 4), 1);
% scatter(repmat(window(1) + 0.5, totalTrials, 1), (1:totalTrials)', 4, c);

axes(hax(2)); hold on;
phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('HDB/ChAT'); set(gca, 'XLim', window, 'YTick', []);

if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(3)); hold on;
    phRasterFromTE(TE, trialsByType{1}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('VTA/DAT'); set(gca, 'XLim', window, 'YTick', []);
end


axes(hax(4)); hold on;
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('HDB/ChAT'); set(gca, 'XLim', window);

if ismember(2, TE.(PhotometryField).settings.channels)
    axes(hax(5)); hold on;
    phRasterFromTE(TE, trialsByType{3}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('VTA/DAT'); set(gca, 'XLim', window, 'YTick', []);
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





    
