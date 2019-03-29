%% averages, aversive, airpuff, also add uncued shock
fdField = 'ZS';
saveName = 'PE_BLA_exampleMouse_aversive_avgs';  
h=ensureFigure(saveName, 1); 
sessionIndexList = 5;


linecolors = [1 0 0; 0 0 0; 1 0 1];            
window = [-7 4];
subplot(2, 1, 1);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Left']);

[ha, hl] = phPlotAverageFromTE(TE, uncuedShock, 1,...
    'zeroTimes', TE.Shock, 'FluorDataField', fdField, 'window', window, 'cmap', [0 0 0]); hold on; %high value, reward