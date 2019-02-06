window = [-4 6];
saveName = 'Aversive_pupilRasters';
ensureFigure(saveName, 1); colormap jet;
subplot(1,3,1);
[CData, XData] = alignedDataWindow(TE, TE.pupil.pupDiameterNorm, true(length(TE.filename), 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Shock, 'window', [-6 4]);
imagesc(CData);