

trialsBack = 20;
cp_licks = bpChangePoints([AR.csMinus.licks_cs.before(goodReversals, end - trialsBack + 1:end) AR.csPlus.licks_cs.after(goodReversals, 1:end)], 2, 1000);

cp_ch1 = bpChangePoints([AR.csMinus.phPeakMean_cs_ch1.before(goodReversals, end - trialsBack + 1:end) AR.csPlus.phPeakMean_cs_ch1.after(goodReversals, 1:end)], 2, 1000);

cp_ch2 = bpChangePoints([AR.csMinus.phPeakMean_cs_ch2.before(goodReversals, end - trialsBack + 1:end) AR.csPlus.phPeakMean_cs_ch2.after(goodReversals, 1:end)], 2, 1000);

ensureFigure('test_cp', 1);
subplot(2,2,1);
scatter(cp_licks.index, cp_licks.logit)
subplot(2,2,2);
scatter(cp_ch1.index, cp_ch1.logit)
subplot(2,2,3);
scatter(cp_ch2.index, cp_ch2.logit)

%% 
criterion = 2;
ensureFigure('test_cp_scatter', 1);
goodrev = (cp_licks.logit > criterion) & (cp_ch1.logit > criterion) & (cp_ch2.logit > criterion);
scatter(cp_ch1.index(goodrev), cp_ch2.index(goodrev)); addUnityLine; 
indexDifference = cp_ch2.index(goodrev) - cp_ch1.index(goodrev);
mean(indexDifference)
[p, h] = signrank(indexDifference)

ensureFigure('cp_cum_dist', 1);
axes; hold on;
[sorted, index] = cum(cp_licks.index(goodrev));
plot(sorted, index, 'k');
[sorted, index] = cum(cp_ch1.index(goodrev));
plot(sorted, index, 'g');
[sorted, index] = cum(cp_ch2.index(goodrev));
plot(sorted, index, 'r');

