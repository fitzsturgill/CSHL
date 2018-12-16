%% plot example reversals with weibull fits, changepoints
nShow = 6;
% good revs include [77 3 48 74];,   3 and 74 best examples
% mediocre,  84 and 91 and 42 are mediocre ones
% 4 is good but gradual
toShow = find(goodReversals);

ordering = randperm(sum(goodReversals), nShow);
toShow = toShow(ordering);
toShow2 = 1:sum(goodReversals);
toShow2 = toShow2(ordering);

% nShow = 4;  toShow = [74 4 91 42];
ensureFigure('reversals_pooled_Latency_examples', 1);

for counter = 1:nShow   
    thisRev = toShow(counter);
    thisRev2 = toShow2(counter);
    % weibull
    subplot(nShow, 3, counter*3 - 2);
%     plot(weibull(thisRev).toFit, 'g.'); hold on;
%     plot(weibull(thisRev).object); legend off;
    set(gca, 'XLim', [0 120]);
    % changepoint
    subplot(nShow, 3, counter*3 - 1); hold on;
    scatter(1:length(cp_licks.cumsum{thisRev}), cp_licks.cumsum{thisRev}, 10, cp_licks.logitAll{thisRev}); colormap jet;
    line(repmat(cp_licks.index(thisRev), 1, 2), get(gca, 'YLim')); set(gca, 'XLim', [0 120]);
    textBox(sprintf('Logit=%.2f', cp_licks.logit(thisRev)));
    % exponential
    subplot(nShow, 3, counter*3); hold on;    
%     plot(expFit.csPlus.licks_cs.toFit{thisRev2}, 'g.'); hold on;
%     plot(expFit.csPlus.licks_cs.object{thisRev2}); legend off;
    set(gca, 'XLim', [0 120]);
end

return;
%%
    ensureFigure('test', 1);
    scatter(1:length(cp_licks.cumsum{thisRev}), cp_licks.cumsum{thisRev}, 10, cp_licks.logitAll{thisRev}); colormap jet;
    line(repmat(cp_licks.index(thisRev), 1, 2), get(gca, 'YLim')); set(gca, 'XLim', [0 120]);
    textBox(sprintf('Logit=%.2f', cp_licks.logit(thisRev)));