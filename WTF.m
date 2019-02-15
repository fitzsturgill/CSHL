
ensureFigure('test', 1);

for counter = 1:6
    subplot(3,2,counter); plot(TE.Photometry.xData', [TE.Photometry.data(1).raw(counter,:)'  TE.Photometry.data(2).raw(counter,:)']);
end