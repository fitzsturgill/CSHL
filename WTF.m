ensureFigure('test', 1);
plot(TE.Photometry.bleachFit(counter, channel).templateX, TE.Photometry.bleachFit(counter, channel).trialTemplate, 'g'); hold on;      
plot(TE.Photometry.bleachFit(counter, channel).fitX, TE.Photometry.bleachFit(counter, channel).trialTemplateFull, 'b'); 