ensureFigure('raw', 1);

% nTrials = size(sessions.SessionData.NidaqData, 1);
trials = [2:8 10:12];
nTrials = length(trials);
data.mean = zeros(nTrials, 1); data.std = zeros(nTrials, 1);
axes; hold on;
for counter = 1:nTrials
    trial = trials(counter);
    subplot(3,4, counter);
    thisdata = SessionData.NidaqData{trial, 1}(6101:54900,2);
    plot(detrend(thisdata));
    set(gca, 'YLim', [-0.02 0.02]);
    data.mean(counter) = mean(thisdata);
    data.std(counter) = std(detrend(thisdata));
    data.data{counter} = thisdata;
end
 
ensureFigure('scatter', 1);
scatter(data.mean, data.std);