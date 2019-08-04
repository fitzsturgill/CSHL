% Ach_BLA_xcorr_wheel_figure

saveOn = 1;
figSize = [2 0.5];
savepath = 'Z:\SummaryAnalyses\ChAT_PE_Manuscript';
Fs = 20;
window = [-20 0];
maxLagInSeconds = 5;
maxLag = round(maxLagInSeconds * Fs);
ylim = [-0.2 0.6];
%% together

saveName = 'ACh_LR_wheel_Xcorr_spont';
ensureFigure(saveName, 1);
loadPath = 'Z:\SummaryAnalyses\wheel_v1\ACh_7_wheel_v1_Feb21_2019_Session1';
load(fullfile(loadPath, 'TE.mat'), 'TE');
[data_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, window);
data_chat = extractDataByTimeStamps(TE.Photometry.data(1).ZS, TE.Photometry.startTime, 20, TE.Reward, window);
[r, shiftR, rawR, lags] = correctedXCorr(data_chat, data_dat, maxLag, 2);

subplot(1,2,2); hold on;
line([0 0], ylim, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
plot(lags * 1/Fs, r, 'k');
% xlabel('Lag (s)');
set(gca, 'YLim', ylim, 'YTick', []);


loadPath = 'Z:\SummaryAnalyses\wheel_v1\ACh_15_wheel_v1_Apr17_2019_Session1';
load(fullfile(loadPath, 'TE.mat'), 'TE');
[data_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, window);
data_chat = extractDataByTimeStamps(TE.Photometry.data(1).ZS, TE.Photometry.startTime, 20, TE.Reward, window);
[r, shiftR, rawR, lags] = correctedXCorr(data_chat, data_dat, maxLag, 2);

subplot(1,2,1); hold on;
line([0 0], ylim, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
plot(lags * 1/Fs, r, 'k'); hold on;
% xlabel('Lag (s)');
set(gca, 'YLim', ylim, 'YTick', [0 0.5]);

formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg'])); 
end