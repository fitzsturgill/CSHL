

%% images
% orderings = {'newCsPlus', 'newCsMinus', 'alwaysCsPlus'};
fieldsToShow = {'phPeakMean_cs_ch1', 'phPeakMean_cs_ch2', 'licks_cs', 'csLicksROC', 'pupil_cs', 'whisk_cs'};
titles = {'ACh', 'Dop.', 'Licks', 'Lick auROC', 'Pupil', 'Whisk'};
clim = [-5 5];
fh=[];

savename = ['newCsPlus_image' expType];
fh(end+1) = ensureFigure(savename, 1);
cLimFactor = 3;
smoothWindow = 1;
xData = [min(newCsPlus_trialNumber), max(newCsPlus_trialNumber)];
xlim = [-30 70];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = newCsPlus.(sfield)(sortOrder, :);
    cData = smoothdata(cData, 2, 'movmean', smoothWindow, 'omitnan');
    imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:length(sortOrder), [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 length(sortOrder)]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs+');
set(gcf, 'Position', [304   217   633   485]);

savename = ['newCsMinus_image' expType];
fh(end+1) = ensureFigure(saveName, 1);
cLimFactor = 3;
xData = [min(newCsMinus_trialNumber), max(newCsMinus_trialNumber)];
xlim = [-30 70];
for fcounter = 1:length(fieldsToShow)
    sfield = fieldsToShow{fcounter};
    subplot(2,3,fcounter);
    cData = newCsMinus.(sfield)(sortOrder, :);
    cData = smoothdata(cData, 2, 'movmean', smoothWindow, 'omitnan');
    imagesc('XData', xData, 'CData', cData); set(gca, 'XLim', xlim); hold on; 
    scatter(zeros(nReversals, 1) + xlim(1) + 1, 1:length(sortOrder), [], repmat(goodReversals(sortOrder), 1, 3) .* [1 0 0], 's', 'filled'); 
    set(gca, 'YLim', [1 length(sortOrder)]);
    set(gca, 'CLim', [nanmean(nanmean(cData, 1), 2) - nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor nanmean(nanmean(cData, 1), 2) + nanstd(nanstd(cData, 0, 1), 0, 2) * cLimFactor]);
    t = textBox(titles{fcounter}, gca, [0.1 0.95]); set(t, 'Color', [0 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end
subplot(2,3,5); xlabel('Odor presentations from reversal');
subplot(2,3,2); title('New Cs-');
set(gcf, 'Position', [304   217   633   485]);



%%

% dT = 0.05;
% 
% sr = 40;
% 
% duration = 60;
% 
% nTrials = 50;
% 
% stimes = cell(nTrials, 1);
% for counter = 1:nTrials
%     stimes{counter} = makeSpikes(dT, sr, duration, 1);
% end
% fitz= 
% binraster = stimes2binraster(stimes, 0:dT:duration, dT, repmat([0 60], nTrials, 1));
% 
% 
% 
% % binraster = stimes2binraster(spikes
% 
% 
% function spikes = makeSpikes(timeStepS, spikesPerS, durationS, numTrains)
% 
% if (nargin < 4)
%     numTrains = 1;
% end
% times = [0:timeStepS:durationS];
% spikes = zeros(numTrains, length(times));
% for train = 1:numTrains
%     vt = rand(size(times));
%     spikes(train, :) = (spikesPerS*timeStepS) > vt;
% end
% end