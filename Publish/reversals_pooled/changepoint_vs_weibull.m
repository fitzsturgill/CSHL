    
    
    
    
    %% cross correlations, new Cs+, corrected

savename = 'xcorr_newCsPlus_corrected_perm';
ensureFigure(savename, 1); 
% [R, lags] = avgXCorr(x, y, maxlag)
maxlag = 10;
trialRange = [30 30];
R = deal(zeros(maxlag*2+1, 6, 3)); % corrected, shift predictor, raw occupy third dimension
lags = -maxlag:maxlag;

%     chat = nanzscore(newCsPlus.phPeakMean_cs_ch1, 0, 2);
%     dat = nanzscore(newCsPlus.phPeakMean_cs_ch2, 0, 2);
%     licks = nanzscore(newCsPlus.licks_cs, 0, 2);

     

chat = perm_ch1(:, newCsPlus.firstRevTrial - trialRange:newCsPlus.firstRevTrial + trialRange - 1);
dat = perm_ch2(:, newCsPlus.firstRevTrial - trialRange:newCsPlus.firstRevTrial + trialRange - 1);
licks = perm_licks(:, newCsPlus.firstRevTrial - trialRange:newCsPlus.firstRevTrial + trialRange - 1);  



[R1, R2, R3, testlags] = correctedXCorr(chat, dat, maxlag, 2);
R(:,1,1) = R1; R(:,1,2) = R2; R(:,1,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(chat, chat, maxlag, 2);
R(:,2,1) = R1; R(:,2,2) = R2; R(:,2,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(dat, dat, maxlag, 2);
R(:,3,1) = R1; R(:,3,2) = R2; R(:,3,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(chat, licks, maxlag, 2);
R(:,4,1) = R1; R(:,4,2) = R2; R(:,4,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(dat, licks, maxlag, 2);
R(:,5,1) = R1; R(:,5,2) = R2; R(:,5,3) = R3;

[R1, R2, R3, ~] = correctedXCorr(licks, licks, maxlag, 2);
R(:,6,1) = R1; R(:,6,2) = R2; R(:,6,3) = R3;
    
   

subplot(3,2,1); plot(lags, squeeze(R(:,1,:))); title('chat vs dat'); legend('corrected', 'shift predictor', 'raw'); legend boxoff;
subplot(3,2,2); plot(lags, squeeze(R(:,2,:))); title('chat');  legend('corrected', 'shift predictor', 'raw'); legend boxoff;
subplot(3,2,3); plot(lags, squeeze(R(:,3,:))); title('dat'); 
subplot(3,2,4); plot(lags, squeeze(R(:,4,:))); title('chat vs licks');
subplot(3,2,5); plot(lags, squeeze(R(:,5,:))); title('dat vs licks'); xlabel('new cs+ trials');
subplot(3,2,6); plot(lags, squeeze(R(:,6,:))); title('licks'); xlabel('new cs+ trials');


if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end



%% align reversals by lick tt20 value computed from Weibull fit

trialWindow = [-20 20];
cpField = 'licks_cs';
cLimFactor = 2;
trialWindow_rev = [0 40];
markerSize = 4;


% first new csPlus (acquisition)

% setup for aligned by changepoint
goodOnes = goodReversals & (weibull.csPlus.(cpField).latency > 0);
zeroTrials = weibull.csPlus.(cpField).latency(goodOnes);
reversalPoints = 0 - zeroTrials;
reversalPoints = reversalPoints;
[sorted, ix] = sort(reversalPoints); % THEN sort them


% first select good subset
good_subtract = newCsPlus.phPeakMean_cs_AchMinusDop(goodOnes, :);
good_licks = newCsPlus.licks_cs(goodOnes, :);
% try normalizing licks
good_licks = good_licks ./ percentile(good_licks(:,newCsPlus.firstRevTrial:end), 0.9, 2);
good_ch1 = newCsPlus.phPeakMean_cs_ch1(goodOnes, :);
good_ch2 = newCsPlus.phPeakMean_cs_ch2(goodOnes, :);

trials = true(sum(goodOnes), 1);
[aligned_subtract, xData] = alignedDataWindow(good_subtract, trials, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1), sum(goodOnes), 1));
[aligned_licks, ~] = alignedDataWindow(good_licks, trials, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch1, ~] = alignedDataWindow(good_ch1, trials, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
[aligned_ch2, ~] = alignedDataWindow(good_ch2, trials, 'zeroTimes', zeroTrials, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));

% randomly permute across reversals (not time) to preserve average but
% scramble potential changepoints
maxTrials = length(newCsPlus.trialNumber);

% generate subscript indices for permutation
row_ix = repmat(1:maxTrials, sum(goodOnes), 1);
col_ix = zeros(sum(goodOnes), maxTrials);
for counter = 1:maxTrials
    col_ix(:,counter) = randperm(sum(goodOnes))';
end

% lin_ix = sub2ind([sum(goodOnes) maxTrials], col_ix, row_ix);
% perm_subtract = smoothdata(good_subtract(lin_ix), 2, 'movmean', smoothWindow, 'omitnan');
% perm_licks = good_licks(lin_ix);
% perm_licks = smoothdata(perm_licks, 2, 'movmean', smoothWindow, 'omitnan');
% perm_ch1 = smoothdata(good_ch1(lin_ix), 2, 'movmean', smoothWindow, 'omitnan');
% perm_ch2 = smoothdata(good_ch2(lin_ix), 2, 'movmean', smoothWindow, 'omitnan');
% 
% % quick plot of regular and permuted lick_cs images, and corresponding
% % averages
% savename = 'permutation_sanity_check';
% ensureFigure(savename, 1);
% clim = [0 1]; xlim = [-20 40];
% subplot(3,1,1);
% imagesc(newCsPlus.trialNumber(:, newCsPlus.firstRevTrial - 21:newCsPlus.firstRevTrial + 39), 1:sum(goodOnes), good_licks(:, newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40), clim);
% set(gca, 'XLim', xlim);
% title('Licks');
% subplot(3,1,2);
% imagesc(newCsPlus.trialNumber(:, newCsPlus.firstRevTrial - 21:newCsPlus.firstRevTrial + 39), 1:sum(goodOnes), perm_licks(:, newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40), clim);
% set(gca, 'XLim', xlim);
% title('Licks permuted');
% subplot(3,1,3); hold on; 
% title('averages');
% plot(newCsPlus.trialNumber(newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40), nanmean(good_licks(:, newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40)), '-k');
% plot(newCsPlus.trialNumber(newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40), nanmean(perm_licks(:, newCsPlus.firstRevTrial - 20:newCsPlus.firstRevTrial + 40)), '--', 'Color', [0.8 0.8 0.8]);
% legend({'intact', 'permuted'}, 'Location', 'best');
% set(gca, 'XLim', xlim); xlabel('new Cs+ trials from rev.');
% 
% formatFigurePoster([3 5], '', 10);
% 
% if saveOn
%     saveas(gcf, fullfile(savepath, [savename '.fig']));
%     saveas(gcf, fullfile(savepath, [savename '.jpg']));   
% end    



% cp_perm.csPlus.licks_cs = bpChangePoints(perm_licks(:,newCsPlus.firstRevTrial - baselineTrials:end), 2, 1000, 'up');
% zeroTrials_perm = cp_perm.csPlus.licks_cs.index - baselineTrials;
% reversalPoints_perm = 0 - zeroTrials_perm;
% reversalPoints_perm = reversalPoints_perm;
% [sorted_perm, ix_perm] = sort(reversalPoints_perm); % THEN sort them
% 
% [aligned_subtract_perm, xData] = alignedDataWindow(perm_subtract, trials, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1), sum(goodOnes), 1));
% [aligned_licks_perm, ~] = alignedDataWindow(perm_licks, trials, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
% [aligned_ch1_perm, ~] = alignedDataWindow(perm_ch1, trials, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
% [aligned_ch2_perm, ~] = alignedDataWindow(perm_ch2, trials, 'zeroTimes', zeroTrials_perm, 'window', trialWindow, 'Fs', 1, 'startTimes', repmat(newCsPlus.trialNumber(1),  sum(goodOnes), 1));
% 
% 
% % setup for alignment by reversal
% % [sorted_rev, sortOrder_rev] = sort(cp.csPlus.licks_cs.index(goodOnes));
cp_rev = zeroTrials;
% cp_rev_perm = zeroTrials_perm;

% first images, sorted by reversal point.
savename = ['tt20_aligned_newCsPlus_images' '_' expType];
ensureFigure(savename, 1);

ckey = pink;
nColors = size(ckey, 1);

% max_logit = max(cp.csPlus.licks_cs.logit(goodOnes & isfinite(cp.csPlus.licks_cs.logit)));
% max_logit = max_logit * 1.2; % make Inf look brighter than others
% color_ix = cp.csPlus.licks_cs.logit(goodOnes);
% color_ix(~isfinite(color_ix)) = max_logit;
% color_ix = ceil(color_ix ./ max_logit .* nColors);
% cp_colors = ckey(color_ix, :);


% color_ix_perm = cp_perm.csPlus.licks_cs.logit;
% color_ix_perm(~isfinite(color_ix)) = max_logit;
% color_ix_perm = ceil(color_ix ./ max_logit .* nColors);
% cp_colors_perm = ckey(color_ix, :);



subplot(4,4,1); hold on;
title('Reversal Aligned');
ylabel('Cue licks');
cData = good_licks;
% use same clim for all images in each row
clim =  [nanmean(cData(isfinite(cData))) - nanstd(cData(isfinite(cData))) * cLimFactor nanmean(cData(isfinite(cData))) + nanstd(cData(isfinite(cData))) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);

subplot(4,4,2); hold on;
title('trials to 20% Aligned');
cData = aligned_licks(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4);
set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

% subplot(4,4,3); hold on;
% title('Rev. aligned, permuted');
% cData = perm_licks;
% % use same clim for all images in each row
% clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
% imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
% scatter(cp_rev_perm, 1:sum(goodOnes), 20, cp_colors_perm, 'filled');
% set(gca, 'YLim', [1 sum(goodOnes)]);
% 
% 
% subplot(4,4,4); hold on;
% title('CP aligned, permuted');
% cData = aligned_licks_perm(ix_perm, :);
% imagesc('XData', xData, 'CData', cData, clim);
% plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,5); hold on;
ylabel('ACh.', 'Color', mycolors('ChAT'));
cData = newCsPlus.phPeakMean_cs_ch1(goodOnes, :);
% use same clim for all images in each row
clim =  [nanmean(cData(isfinite(cData))) - nanstd(cData(isfinite(cData))) * cLimFactor nanmean(cData(isfinite(cData))) + nanstd(cData(isfinite(cData))) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

subplot(4,4,6); hold on;
cData = aligned_ch1(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

% subplot(4,4,7); hold on;
% cData = perm_ch1;
% % use same clim for all images in each row
% clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
% imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
% scatter(cp_rev_perm, 1:sum(goodOnes), 20, cp_colors_perm, 'filled');
% set(gca, 'YLim', [1 sum(goodOnes)]);
% 
% subplot(4,4,8); hold on;
% cData = aligned_ch1_perm(ix_perm, :);
% imagesc('XData', xData, 'CData', cData, clim);
% plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,9); hold on;
ylabel('Dop.', 'Color', mycolors('DAT'));
cData = newCsPlus.phPeakMean_cs_ch2(goodOnes, :);
% use same clim for all images in each row
clim =  [nanmean(cData(isfinite(cData))) - nanstd(cData(isfinite(cData))) * cLimFactor nanmean(cData(isfinite(cData))) + nanstd(cData(isfinite(cData))) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)], 'XTick', []);

subplot(4,4,10); hold on;
cData = aligned_ch2(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

% subplot(4,4,11); hold on;
% cData = perm_ch2;
% % use same clim for all images in each row
% clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
% imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
% scatter(cp_rev_perm, 1:sum(goodOnes), 20, cp_colors_perm, 'filled');
% set(gca, 'YLim', [1 sum(goodOnes)]);
% 
% subplot(4,4,12); hold on;
% cData = aligned_ch2_perm(ix_perm, :);
% imagesc('XData', xData, 'CData', cData, clim);
% plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', [], 'XTick', []);

subplot(4,4,13); hold on;
ylabel('Subtract', 'Color', [0 1 0]);
cData = newCsPlus.phPeakMean_cs_AchMinusDop(goodOnes, :);
% use same clim for all images in each row
clim =  [nanmean(cData(isfinite(cData))) - nanstd(cData(isfinite(cData))) * cLimFactor nanmean(cData(isfinite(cData))) + nanstd(cData(isfinite(cData))) * cLimFactor];
imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
scatter(cp_rev, 1:sum(goodOnes), 20, [1 1 1], 'filled');
set(gca, 'YLim', [1 sum(goodOnes)]);
xlabel('Trials from reversal');

subplot(4,4,14); hold on;
cData = aligned_subtract(ix, :);
imagesc('XData', xData, 'CData', cData, clim);
plot(reversalPoints(ix), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', []);
xlabel('Trials from tt20');

% subplot(4,4,15); hold on;
% cData = perm_subtract;
% % use same clim for all images in each row
% clim =  [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor];
% imagesc('XData', newCsPlus.trialNumber, 'CData', cData, clim); set(gca, 'XLim', trialWindow_rev); hold on; 
% scatter(cp_rev_perm, 1:sum(goodOnes), 20, cp_colors_perm, 'filled');
% set(gca, 'YLim', [1 sum(goodOnes)]);
% 
% subplot(4,4,16); hold on;
% cData = aligned_subtract_perm(ix_perm, :);
% imagesc('XData', xData, 'CData', cData, clim);
% plot(reversalPoints_perm(ix_perm), 1:sum(goodOnes), ':w', 'LineWidth', 4); set(gca, 'XLim', trialWindow, 'YLim', [1 sum(goodOnes)], 'YTick', []);
% xlabel('Trials from changepoint');



formatFigurePoster([8 10], '', 10);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

% second averages
savename = ['tt20_aligned_newCsPlus_avgs' '_' expType];
ensureFigure(savename, 1);
marker = '.';
subplot(2,2,1); title('new Cs+ (acquisition)');
[hl, hp] = boundedline(xData, [nanmean(aligned_licks)' ], permute([nanSEM(aligned_licks)'], [1 3 2]),...
    'cmap', [0 0 0], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
set(hl, 'Marker', marker);
ylabel('cue licks');
addOrginLines;
legend(hl, 'Intact', 'Permuted', 'Location', 'best'); legend boxoff;

subplot(2,2,2);
[hl, hp] = boundedline(xData, [nanmean(aligned_ch1)' ], permute([nanSEM(aligned_ch1)' ], [1 3 2]),...
    'cmap', [mycolors('ChAT')], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
set(hl, 'Marker', marker);
ylabel('ACh.');
addOrginLines;

subplot(2,2,3);
[hl, hp] = boundedline(xData, [nanmean(aligned_ch2)' ], permute([nanSEM(aligned_ch2)' ], [1 3 2]),...
    'cmap', [mycolors('DAT')], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
set(hl, 'Marker', marker);
ylabel('Dop.');
xlabel('trials from lick tt20');
addOrginLines;

subplot(2,2,4);
[hl, hp] = boundedline(xData, [nanmean(aligned_subtract)' ], permute([nanSEM(aligned_subtract)' ], [1 3 2]),...
    'cmap', [0 1 0;], 'nan', 'gap'); hold on
set(gca, 'XLim', trialWindow);%, 'YLim', [-1 2]);
ylabel('ACh. - Dop.');
xlabel('trials from lick tt20');
addOrginLines;

formatFigurePoster([5.5 4], '', 8);

if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    


