window = [-20 0];
[data_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, window);
data_chat = extractDataByTimeStamps(TE.Photometry.data(1).ZS, TE.Photometry.startTime, 20, TE.Reward, window);


maxLagInSeconds = 5;
maxLag = round(maxLagInSeconds * Fs);

[r, shiftR, rawR, lags] = correctedXCorr(data_chat, data_dat, maxLag, 2);
ensureFigure('xcorr', 1);
plot(lags * 1/Fs, r, 'k'); hold on;
plot(lags * 1/Fs, rawR, 'r');
plot(lags * 1/Fs, shiftR, 'b');
xlabel('Time (s)'); ylabel('ChAT x DAT XCorr'); 
legend({'corrected', 'raw', 'shift predictor'});


if saveOn
    saveas(gcf, fullfile(savepath, 'xcorr_preReward.fig'));
    saveas(gcf, fullfile(savepath, 'xcorr_preReward.jpg'));
end