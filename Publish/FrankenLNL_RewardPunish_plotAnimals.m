
xwindow = [-3 3];
ensureFigure('Reward_vs_LickLatency', 1);

subplot(1,6,1); eventRasterFromTE(TE, rewardTrials & ~uncuedTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroTimes', TE.Outcome, 'window', xwindow, 'sortValues', TE.lickLatency_us);
set(gca, 'XLim', xwindow);
subplot(1,6,2); phRasterFromTE(TE, rewardTrials & ~uncuedTrials, 1, 'FluorDataField', 'ZS', 'sortValues', TE.lickLatency_us	, 'zeroTimes', TE.Us, 'window', [-3 3]);
subplot(1,6,3); phRasterFromTE(TE, rewardTrials & ~uncuedTrials, 2, 'FluorDataField', 'ZS', 'sortValues', TE.lickLatency_us, 'zeroTimes', TE.Us, 'window', [-3 3]);
subplot(1,6,4); eventRasterFromTE(TE, rewardTrials & uncuedTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroTimes', TE.Outcome, 'window', xwindow, 'sortValues', TE.lickLatency_us);
set(gca, 'XLim', xwindow);
subplot(1,6,5); phRasterFromTE(TE, rewardTrials & uncuedTrials, 1, 'FluorDataField', 'ZS', 'sortValues', TE.lickLatency_us, 'zeroTimes', TE.Us, 'window', [-3 3]);
subplot(1,6,6); phRasterFromTE(TE, rewardTrials & uncuedTrials, 2, 'FluorDataField', 'ZS', 'sortValues', TE.lickLatency_us, 'zeroTimes', TE.Us, 'window', [-3 3]);


ensureFigure('Reward_vs_LickRate', 1);

% lickrates = sort(TE.licks_us.rate(rewardTrials & ~uncuedTrials));
% subplot(1,6,1); plot(lickrates, 1:sum(rewardTrials & ~uncuedTrials), '-k'); set(gca, 'YDir', 'reverse');
subplot(1,6,1); eventRasterFromTE(TE, rewardTrials & ~uncuedTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroTimes', TE.Outcome, 'window', xwindow, 'sortValues', TE.licks_us.rate);
set(gca, 'XLim', xwindow);
subplot(1,6,2); phRasterFromTE(TE, rewardTrials & ~uncuedTrials, 1, 'FluorDataField', 'ZS', 'sortValues', TE.licks_us.rate	, 'zeroTimes', TE.Us, 'window', [-3 3]);
subplot(1,6,3); phRasterFromTE(TE, rewardTrials & ~uncuedTrials, 2, 'FluorDataField', 'ZS', 'sortValues', TE.licks_us.rate, 'zeroTimes', TE.Us, 'window', [-3 3]);
lickrates = sort(TE.licks_us.rate(rewardTrials & uncuedTrials));
% subplot(1,6,4); plot(lickrates, 1:sum(rewardTrials & uncuedTrials), '-k'); set(gca, 'YDir', 'reverse');
% set(gca, 'XLim', xwindow);
subplot(1,6,4); eventRasterFromTE(TE, rewardTrials & uncuedTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroTimes', TE.Outcome, 'window', xwindow, 'sortValues', TE.licks_us.rate);
set(gca, 'XLim', xwindow);
subplot(1,6,5); phRasterFromTE(TE, rewardTrials & uncuedTrials, 1, 'FluorDataField', 'ZS', 'sortValues', TE.licks_us.rate, 'zeroTimes', TE.Us, 'window', [-3 3]);
subplot(1,6,6); phRasterFromTE(TE, rewardTrials & uncuedTrials, 2, 'FluorDataField', 'ZS', 'sortValues', TE.licks_us.rate, 'zeroTimes', TE.Us, 'window', [-3 3]);