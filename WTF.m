%% signal correlations

xData = [us_pooled.rew.avg_delta(:,1) us_pooled.puff.avg_delta(:,1) us_pooled.shock.avg_delta(:,1)];
yData = [us_pooled.rew.avg_delta(:,2) us_pooled.puff.avg_delta(:,2) us_pooled.shock.avg_delta(:,2)];

Rsignal = zeros(size(xData, 1));

for counter = 1:size(xData, 1)
    Rsignal(counter) = corr(xData(counter, :), yData(counter, :));
end