PhotometryField = 'Photometry';
trialNumbering = 'global';
CLimFactor = 3;
window = [-4 6];

saveName = 'cued_reinforcment_allBehavior_shock_special';
ensureFigure(saveName, 1);

subplot(1,8,1);
eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
    'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('cued reward');         
addStimulusPatch(gca, [0 1]);

subplot(1,8,2);
phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('HDB'); set(gca, 'XLim', window);

subplot(1,8,3);
phRasterFromTE(TE, trialsByType{1}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('HDB'); set(gca, 'XLim', window);



subplot(1,8,4);
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('Cued Shock: HDB/ChAT'); set(gca, 'XLim', window);

subplot(1,8,5);
phRasterFromTE(TE, trialsByType{3}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('VTA/DAT'); set(gca, 'XLim', window);

subplot(1,8,6);
alignedDataRaster(TE.pupil.pupDiameterNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate);
title('Pupil diameter'); set(gca, 'XLim', window);

subplot(1,8,7);
alignedDataRaster(TE.pupil.frameAvgNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate);
title('Eye closure'); set(gca, 'XLim', window);

subplot(1,8,8);
alignedDataRaster(TE.Wheel.data.V, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.Wheel.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.Wheel.Fs);
title('Eye closure'); set(gca, 'XLim', window);