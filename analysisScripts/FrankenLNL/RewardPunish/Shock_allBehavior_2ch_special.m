saveName = 'cued_reinforcement_allBehavior_shock_special';
fig = ensureFigure(saveName, 1);
mcLandscapeFigSetup(fig);

%% Define the positions of different axes matrices on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    matpos_title_1 = [0 0 3/8 .1];
    matpos_title_2 = [3/8 0 (1 - 3/8) .1];
    matpos_rasters = [0 .1 1 1];    
    params.cellmargin = [.05 .05 0.05 0.05];    
    
%% add titles
params.matpos = matpos_title_1;
ax = textAxes(fig, 'Cued Reward', 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
params.matpos = matpos_title_2;
ax = textAxes(fig, 'Cued Shock', 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
    

%% cued Reward, cued Shock, all behavior
PhotometryField = 'Photometry';
trialNumbering = 'global';
CLimFactor = 3;
window = [-4 6];
totalTrials = length(TE.filename);


params.matpos = matpos_rasters;


ax = axes; hold on;
eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
    'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('Cued Reward:');         
addStimulusPatch(gca, [0 2]);

% add scatter plot with tags for block type blue = 
c = zeros(totalTrials, 3);
c(TE.BlockNumber == 2, :) = repmat([0 0 1], sum(TE.BlockNumber == 2), 1);
c(TE.BlockNumber == 6, :) = repmat([0.5 0.5 0.5], sum(TE.BlockNumber == 6), 1);
c(TE.BlockNumber == 4, :) = repmat([1 0 0], sum(TE.BlockNumber == 4), 1);
scatter(repmat(window(1) + 0.5, totalTrials, 1), (1:totalTrials)', 4, c);

ax(2) = subplot(1,8,2);
phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('HDB/ChAT'); set(gca, 'XLim', window, 'YTick', []);


ax(3) = subplot(1,8,3);
phRasterFromTE(TE, trialsByType{1}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('VTA/DAT'); set(gca, 'XLim', window, 'YTick', []);



ax(4) = subplot(1,8,4);
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('Cued Shock: HDB/ChAT'); set(gca, 'XLim', window);

ax(5) = subplot(1,8,5);
phRasterFromTE(TE, trialsByType{3}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('VTA/DAT'); set(gca, 'XLim', window, 'YTick', []);

ax(6) = subplot(1,8,6);
alignedDataRaster(TE.pupil.pupDiameterNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Pupil diameter'); set(gca, 'XLim', window, 'YTick', []);

ax(7) = subplot(1,8,7);
alignedDataRaster(TE.pupil.frameAvgNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Eye closure'); set(gca, 'XLim', window, 'YTick', []);

ax(8) = subplot(1,8,8);
alignedDataRaster(TE.Wheel.data.V, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.Wheel.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.Wheel.Fs, 'CLim', [0 20]);
title('Velocity'); set(gca, 'XLim', window, 'YTick', []);

for counter = 1:length(ax)
    line(ax(counter), [0 0], get(ax(counter), 'YLim'), 'Color', 'r'); line(ax(counter), [2 2], get(ax(counter), 'YLim'), 'Color', 'r');
end
if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end





    
