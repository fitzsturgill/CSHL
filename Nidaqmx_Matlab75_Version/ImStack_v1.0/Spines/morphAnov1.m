function out = morphAnov1(data,group,type)

%[p2d,tbl2d,stat2d] = anova1(data,group);
[p2d,tbl2d,stat2d] = kruskalwallis(data,group);
figure('pos', [232   130   661   548]);
[c2d,m2d] = multcompare(stat2d);
pairs = showSignificantPairsN(real(c2d),type);
title(['Analysis of ' type ' ' ]);
setLineColors(gca,'black');
copyToClip(gcf);
out = p2d;