ensureFigure('test2', 1);
xlim = [40 70] - 20;
trial = 5;
plot(TE.Photometry.xData, [TE.Photometry.data(1).ZS(trial,:)' TE.Photometry.data(2).ZS(trial,:)']);

set(gca, 'XLim', xlim);