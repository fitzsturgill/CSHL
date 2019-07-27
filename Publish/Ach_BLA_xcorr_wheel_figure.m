% Ach_BLA_xcorr_wheel_figure


loadPath = 'Z:\SummaryAnalyses\wheel_v1\ACh_7_wheel_v1_Feb21_2019_Session1';
savepath = 'Z:\SummaryAnalyses\ChAT_PE_Manuscript';

load(fullfile(loadPath, 'TE.mat'), 'TE');
Fs = 20;
window = [-20 0];
[data_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, window);
data_chat = extractDataByTimeStamps(TE.Photometry.data(1).ZS, TE.Photometry.startTime, 20, TE.Reward, window);


maxLagInSeconds = 5;
maxLag = round(maxLagInSeconds * Fs);

[r, shiftR, rawR, lags] = correctedXCorr(data_chat, data_dat, maxLag, 2);
saveName = 'ACh_7_wheel_Xcorr_spont';
ensureFigure(saveName, 1);
plot(lags * 1/Fs, r, 'k'); hold on;
% plot(lags * 1/Fs, rawR, 'r');
% plot(lags * 1/Fs, shiftR, 'b');
xlabel('Time (s)'); ylabel('R'); 
% legend({'corrected', 'raw', 'shift predictor'});

formatFigurePublish('size', [1.6 1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

