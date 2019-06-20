%% averages


trialWindow = [-20 40];
xlim = trialWindow;
ylim = [-1 2];

% new cs plus
common = (newCsPlus.firstRevTrial + trialWindow(1)):(newCsPlus.firstRevTrial + trialWindow(2));
savename = 'reversals_newCsPlus_zscored';

ch1Data = nanzscore(newCsPlus.phPeakMean_cs_ch1(goodReversals, common), 0, 2);
ch2Data = nanzscore(newCsPlus.phPeakMean_cs_ch2(goodReversals, common), 0, 2);
lickData = nanzscore(newCsPlus.licks_cs(goodReversals, common), 0, 2);

% fit exponentials to the data
fo = fitoptions('Method', 'NonlinearLeastSquares',...
    'Upper', [Inf  Inf Inf],...
    'Lower', [-Inf -Inf 0],...    
    'StartPoint', [0 1 10]...
    );

expFitMeans.csPlus = struct(); 
expModel = 'a + b * exp(-x/c)';
ft = fittype(expModel, 'options', fo);
startFit = bpX2pnt(0, 1, trialWindow(1));
xData = 0:(trialWindow(2));
taus = zeros(1,3);
taus_int = zeros(2,3);
[fitobject, gof, output] = fit(xData', nanmean(ch1Data(:,startFit:end))', ft, fo);
expFitMeans.csPlus.ch1 = fitobject;
taus(1) = fitobject.c;
c = confint(fitobject); taus_int(:,1) = c(:,3);

[fitobject, gof, output] = fit(xData', nanmean(ch2Data(:,startFit:end))', ft, fo);
expFitMeans.csPlus.ch2 = fitobject;
taus(2) = fitobject.c;
c = confint(fitobject); taus_int(:,2) = c(:,3);

[fitobject, gof, output] = fit(xData', nanmean(lickData(:,startFit:end))', ft, fo);
expFitMeans.csPlus.licks = fitobject;
taus(3) = fitobject.c;
c = confint(fitobject); taus_int(:,3) = c(:,3);

% stupid bar graph
ensureFigure('test', 1);
errorbar(1:3, taus, taus_int(1,:), taus_int(2,:));
set(gca, 'XLim', [-1 4]);

fh(end + 1) = ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(ch2Data), nanSEM(ch1Data)',...
    'cmap', mycolors('dat'), 'nan', 'gap'); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(ch1Data), nanSEM(ch1Data)',...
    'cmap', mycolors('chat'), 'nan', 'gap');
hla(2) = hl;
[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(lickData), nanSEM(lickData)',...
    'cmap', mycolors('licks'), 'nan', 'gap');
hla(3) = hl;
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
    },...
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');

title('New Cs+');
xlabel('Odor presentations from reversal');
ylabel('Cue response (z-scored)');
formatFigurePoster([5.5 4], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    
%%
ensureFigure('test', 1);
hl = plot(expFitMeans.csPlus.ch1, 'predfunc'); hold on;
set(hl, 'Color', mycolors('chat'));
hl = plot(expFitMeans.csPlus.ch2, 'predfunc');
set(hl, 'Color', mycolors('dat'));
plot(expFitMeans.csPlus.licks, 'predfunc');
set(hl, 'Color', mycolors('licks'));
%%

% new cs minus
common = (newCsMinus.firstRevTrial + trialWindow(1)):(newCsMinus.firstRevTrial + trialWindow(2) - 1);
savename = 'reversals_newCsMinus_zscored';
fh(end + 1) = ensureFigure(savename, 1);
hla = zeros(1,3);
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(nanzscore(newCsMinus.phPeakMean_cs_ch2(goodReversals, common), 0, 2)), nanSEM(nanzscore(newCsMinus.phPeakMean_cs_ch2(goodReversals, common), 0, 2))',...
    'cmap', [237 125 49]/256, 'nan', 'gap'); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(nanzscore(newCsMinus.phPeakMean_cs_ch1(goodReversals, common), 0, 2)), nanSEM(nanzscore(newCsMinus.phPeakMean_cs_ch1(goodReversals, common), 0, 2))',...
    'cmap', [171 55 214]/256, 'nan', 'gap');
hla(2) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(nanzscore(newCsMinus.licks_cs(goodReversals, common), 0, 2)), nanSEM(nanzscore(newCsMinus.licks_cs(goodReversals, common), 0, 2))',...
    'cmap', [0.5 0.5 0.5]/256, 'nan', 'gap');
hla(3) = hl;
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
    },...
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');

title('New Cs-');
xlabel('Odor presentations from reversal');
ylabel('Cue response (z-scored)');
formatFigurePoster([5.5 4], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end

%
savename = 'reversals_newCsPlus_subtract_zscored';
ensureFigure(savename, 1);

subtractData = nanzscore(newCsPlus.phPeakMean_cs_ch1(goodReversals, common), 0, 2) - nanzscore(newCsPlus.phPeakMean_cs_ch2(goodReversals, common), 0, 2);

[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(subtractData), nanSEM(subtractData)',...
    'cmap', [0 1 0]); hold on
set(hl, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);    
h  = addOrginLines;
set(hl, 'LineWidth', 2);
legend(hl, 'Ach. - Dop.'); 
title('New Cs+');
xlabel('Odor presentations from reversal');
ylabel('delta ZScore');    
formatFigurePoster([5.5 4], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    


savename = 'reversals_newCsMinus_subtract_zscored';
ensureFigure(savename, 1);

subtractData = nanzscore(newCsMinus.phPeakMean_cs_ch1(goodReversals, common), 0, 2) - nanzscore(newCsMinus.phPeakMean_cs_ch2(goodReversals, common), 0, 2);

[hl, hp] = boundedline(newCsPlus_trialNumber(common), nanmean(subtractData), nanSEM(subtractData)',...
    'cmap', [0 1 0]); hold on
set(hl, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);    
h  = addOrginLines;
set(hl, 'LineWidth', 2);
legend(hl, 'Ach. - Dop.'); 
title('New Cs-');
xlabel('Odor presentations from reversal');
ylabel('delta ZScore');    
formatFigurePoster([5.5 4], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

%%
% new cs minus
common = sum(~isnan(newCsMinus.licks_cs)) > 3;

savename = 'reversals_newCsMinus';
fh(end + 1) = ensureFigure(savename, 1);
hla = zeros(1,6);
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_ch2(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_ch2(goodReversals, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.phPeakMean_cs_ch1(goodReversals, common)), nanSEM(newCsMinus.phPeakMean_cs_ch1(goodReversals, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(newCsMinus_trialNumber(common), nanmean(newCsMinus.licks_cs(goodReversals, common)), nanSEM(newCsMinus.licks_cs(goodReversals, common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch2(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch1(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.licks_cs(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);    
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
            '\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\color[rgb]{0.5,0.5,0.5}Odor 3'},...
            'Location', 'southwest', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
title('New Cs-');
xlabel('Odor presentations from reversal');
ylabel('Cue response');    
formatFigurePoster([5.5 4], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

% always cs plus
common = sum(~isnan(alwaysCsPlus.licks_cs)) > 3;
savename = 'reversals_alwaysCsPlus';
fh(end + 1) = ensureFigure(savename, 1);
hla = zeros(1,6);
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus.phPeakMean_cs_ch2(goodReversals, common)), nanSEM(alwaysCsPlus.phPeakMean_cs_ch2(goodReversals, common))',...
    'cmap', [237 125 49]/256); hold on
hla(1) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus.phPeakMean_cs_ch1(goodReversals, common)), nanSEM(alwaysCsPlus.phPeakMean_cs_ch1(goodReversals, common))',...
    'cmap', [171 55 214]/256);
hla(2) = hl;
[hl, hp] = boundedline(alwaysCsPlus_trialNumber(common), nanmean(alwaysCsPlus.licks_cs(goodReversals, common)), nanSEM(alwaysCsPlus.licks_cs(goodReversals, common))',...
    'cmap', [0.5 0.5 0.5]/256);
hla(3) = hl;
hla(4) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch2(goodReversals, common_odor3)), '--', 'Color', [237 125 49]/256);
hla(5) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.phPeakMean_cs_ch1(goodReversals, common_odor3)), '--', 'Color', [171 55 214]/256);
hla(6) = plot(odor3_trialNumber(common_odor3), nanmean(odor3.licks_cs(goodReversals, common_odor3)), '--', 'Color', [0.5 0.5 0.5]/256);
set(hla, 'LineWidth', 2);
set(gca, 'XLim', xlim);%, 'YLim', ylim);    
h  = addOrginLines;
set(h, 'LineWidth', 2);
legend(hla, {'\bf\color[rgb]{0.9258,0.4883,0.1914}Dop.', '\bf\color[rgb]{0.6680,0.2148,0.8359}Ach.', '\bf\color[rgb]{0.5,0.5,0.5}Licks',...
    '\color[rgb]{0.9258,0.4883,0.1914}Odor 3', '\color[rgb]{0.6680,0.2148,0.8359}Odor 3', '\color[rgb]{0.5,0.5,0.5}Odor 3'},...
    'Location', 'best', 'FontSize', 12, 'Interpreter', 'tex', 'Box', 'off');
title('Always Cs+');
xlabel('Odor presentations from reversal');
ylabel('Cue response');    
formatFigurePoster([5.5 4], '', 12);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end