function barErrorPlot(means,sems,labels,graphTitle)

figure();clf; hold on;

figNum = bar(means,'FaceColor',[0.5 0.5 0.5],...
   'EdgeColor',[0.1,0.1,0.1])

errorbar(1:length(means),means,sems,...
    '.','Color',[0.1,0.1,0.1],'LineWidth',1.2,'MarkerSize',0.1);

set(gca(),'XTick',1:length(means))
set(gca(),'XTickLabel',labels,'fontsize',12,'fontweight','b')
title(graphTitle,'fontsize',12,'fontweight','b')