%% cross correlation
duration = 30;
Fs = 20;
trials = 30;

pinkA = zeros(Fs * duration, trials);
pinkB = zeros(Fs * duration, trials);
for counter = 1:trials
    pinkA(:,counter) = pinknoise(duration * Fs);
    pinkB(:,counter) = pinknoise(duration * Fs);
end


maxLagInSeconds = 10;
Fs = 20;
maxLag = round(maxLagInSeconds * Fs);
% [r, lags] = xcorr(data_chat(:,trial), data_dat(:,trial), maxLag);
[r, shiftR, rawR, lags] = correctedXCorr(pinkA, pinkB, maxLag);
ensureFigure('xcorr', 1);
subplot(2,2,1);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('pink XCorr'); 
legend({'corrected', 'raw', 'shift predictor'});

[r, shiftR, rawR, lags] = correctedXCorr(pinkA, pinkA, maxLag);
subplot(2,2,2);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('pink autoCorr'); 
legend({'corrected', 'raw', 'shift predictor'});

whiteA = rand(Fs * duration, trials) - 0.5;
whiteB = rand(Fs * duration, trials) - 0.5;

[r, shiftR, rawR, lags] = correctedXCorr(whiteA, whiteB, maxLag);
subplot(2,2,3);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('white XCorr'); 
legend({'corrected', 'raw', 'shift predictor'});

[r, shiftR, rawR, lags] = correctedXCorr(whiteA, whiteA, maxLag);
subplot(2,2,4);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('white autoCorr'); 
legend({'corrected', 'raw', 'shift predictor'});