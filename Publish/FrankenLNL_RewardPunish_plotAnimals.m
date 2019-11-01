% FrankenLNL_RewardPunish_plotAnimals


DB = dbLoadExperiment('FrankenLNL_RewardPunish');

PhotometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
trialNumbering = 'consecutive';
CLimFactor = 3;
fhr = [];
fha = [];


% load average data
savepath = fullfile(DB.path, ['pooled' filesep]);
load(fullfile(savepath, 'grandAverages.mat'), 'gAvg');
disp(['*** loading: ' fullfile(savepath, 'grandAverages.mat') ' ***']);
load(fullfile(savepath, 'grandAveragesNorm.mat'), 'gAvgNorm');
disp(['*** loading: ' fullfile(savepath, 'grandAveragesNorm.mat') ' ***']);

xlim = [-6 4];
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    
    
    success = dbLoadAnimal(DB, animal);
    display(animal);    
    savepath = fullfile(DB.path, 'animals', animal, filesep);
    
    saveName = sprintf('Averages_%s', animal);  
    fha(end + 1) =ensureFigure(saveName, 1);
    mcPortraitFigSetup(gcf);
    ah = zeros(3,2);
    ah(1,1) = subplot(3,2,1); title('Left BLA'); ylabel('Reward'); textBox(animal); hold on;
    xData = [gAvgNorm.phCue.cuedReward.xData gAvgNorm.phUs.cuedReward.xData];
    
    yData = [gAvgNorm.phCue.omitReward.data(counter * 2 - 1, :) gAvgNorm.phUs.omitReward.data(counter * 2 - 1, :);...
        gAvgNorm.phCue.uncuedReward.data(counter * 2 - 1, :) gAvgNorm.phUs.uncuedReward.data(counter * 2 - 1, :);
        gAvgNorm.phCue.cuedReward.data(counter * 2 - 1, :) gAvgNorm.phUs.cuedReward.data(counter * 2 - 1, :)];
    plot(xData', yData'); set(gca, 'XLim', xlim);
    legend({'omit', 'uncued', 'cued'}); legend boxoff;
    
    ah(1,2) = subplot(3,2,2); title('Right BLA'); hold on;
    xData = [gAvgNorm.phCue.cuedReward.xData gAvgNorm.phUs.cuedReward.xData];
    
    yData = [gAvgNorm.phCue.omitReward.data(counter * 2, :) gAvgNorm.phUs.omitReward.data(counter * 2, :);...
        gAvgNorm.phCue.uncuedReward.data(counter * 2, :) gAvgNorm.phUs.uncuedReward.data(counter * 2, :);
        gAvgNorm.phCue.cuedReward.data(counter * 2, :) gAvgNorm.phUs.cuedReward.data(counter * 2, :)];    
    
    plot(xData', yData'); set(gca, 'XLim', xlim);
    
    % air puff
    ah(2,1) = subplot(3,2,3); ylabel('Air puff'); hold on;
    xData = [gAvgNorm.phCue.cuedPuff.xData gAvgNorm.phUs.cuedPuff.xData];
    
    yData = [gAvgNorm.phCue.omitPuff.data(counter * 2 - 1, :) gAvgNorm.phUs.omitPuff.data(counter * 2 - 1, :);...
        gAvgNorm.phCue.uncuedPuff.data(counter * 2 - 1, :) gAvgNorm.phUs.uncuedPuff.data(counter * 2 - 1, :);
        gAvgNorm.phCue.cuedPuff.data(counter * 2 - 1, :) gAvgNorm.phUs.cuedPuff.data(counter * 2 - 1, :)];
    plot(xData', yData'); set(gca, 'XLim', xlim);
    
    ah(2,2) = subplot(3,2,4);  hold on;
    xData = [gAvgNorm.phCue.cuedPuff.xData gAvgNorm.phUs.cuedPuff.xData];
    
    yData = [gAvgNorm.phCue.omitPuff.data(counter * 2, :) gAvgNorm.phUs.omitPuff.data(counter * 2, :);...
        gAvgNorm.phCue.uncuedPuff.data(counter * 2, :) gAvgNorm.phUs.uncuedPuff.data(counter * 2, :);
        gAvgNorm.phCue.cuedPuff.data(counter * 2, :) gAvgNorm.phUs.cuedPuff.data(counter * 2, :)];    
    
    plot(xData', yData'); set(gca, 'XLim', xlim);
    
    % shock
    ah(3,1) = subplot(3,2,5); ylabel('Shock'); hold on; xlabel('time from Us (s)');
    xData = [gAvgNorm.phCue.cuedShock.xData gAvgNorm.phUs.cuedShock.xData];
    
    yData = [gAvgNorm.phCue.omitShock.data(counter * 2 - 1, :) gAvgNorm.phUs.omitShock.data(counter * 2 - 1, :);...
        gAvgNorm.phCue.uncuedShock.data(counter * 2 - 1, :) gAvgNorm.phUs.uncuedShock.data(counter * 2 - 1, :);
        gAvgNorm.phCue.cuedShock.data(counter * 2 - 1, :) gAvgNorm.phUs.cuedShock.data(counter * 2 - 1, :)];
    plot(xData', yData'); set(gca, 'XLim', xlim);
    
    ah(3,2) = subplot(3,2,6); xlabel('time from Us (s)'); hold on;
    xData = [gAvgNorm.phCue.cuedShock.xData gAvgNorm.phUs.cuedShock.xData];
    
    yData = [gAvgNorm.phCue.omitShock.data(counter * 2, :) gAvgNorm.phUs.omitShock.data(counter * 2, :);...
        gAvgNorm.phCue.uncuedShock.data(counter * 2, :) gAvgNorm.phUs.uncuedShock.data(counter * 2, :);
        gAvgNorm.phCue.cuedShock.data(counter * 2, :) gAvgNorm.phUs.cuedShock.data(counter * 2, :)];    
    
    plot(xData', yData'); set(gca, 'XLim', xlim);
    sameYScale(ah(1,:)); sameYScale(ah(2,:)); sameYScale(ah(3,:));
        

end




    
    
h = waitbar(0, 'slowly writing avg pdfs');pdfavg = fullfile(DB.path, 'pooled', 'FrankenLNL_RewardPunish_avgsNorm_allAnimals.pdf');
for counter = 1:length(fha)    
    if counter == 1
        export_fig(fha(counter),pdfavg);  % write to pdf

    else
        export_fig(fha(counter),'-append',pdfavg);  % write to pdf
    end
    waitbar(counter/(length(fha)));
end
close(h);

% h = waitbar(0, 'slowly writing raster pdfs');
% pdfavg = fullfile(DB.path, 'pooled', 'FrankenLNL_varyRewardSize_rasters_allAnimals.pdf');
% for counter = 1:length(fhr)    
%     if counter == 1
%         export_fig(fhr(counter),pdfavg);  % write to pdf
%     else
%         export_fig(fhr(counter),'-append',pdfavg);  % write to pdf
%     end
%     waitbar(counter/(length(fhr)));
% end
% close(h);

%%
% 
% xwindow = [-3 3];
% ensureFigure('Reward_vs_LickLatency', 1);
% 
% subplot(1,6,1); eventRasterFromTE(TE, rewardTrials & ~uncuedTrials, 'Port1In', 'trialNumbering', 'consecutive',...
%     'zeroTimes', TE.Outcome, 'window', xwindow, 'sortValues', TE.lickLatency_us);
% set(gca, 'XLim', xwindow);
% subplot(1,6,2); phRasterFromTE(TE, rewardTrials & ~uncuedTrials, 1, 'FluorDataField', 'ZS', 'sortValues', TE.lickLatency_us	, 'zeroTimes', TE.Us, 'window', [-3 3]);
% subplot(1,6,3); phRasterFromTE(TE, rewardTrials & ~uncuedTrials, 2, 'FluorDataField', 'ZS', 'sortValues', TE.lickLatency_us, 'zeroTimes', TE.Us, 'window', [-3 3]);
% subplot(1,6,4); eventRasterFromTE(TE, rewardTrials & uncuedTrials, 'Port1In', 'trialNumbering', 'consecutive',...
%     'zeroTimes', TE.Outcome, 'window', xwindow, 'sortValues', TE.lickLatency_us);
% set(gca, 'XLim', xwindow);
% subplot(1,6,5); phRasterFromTE(TE, rewardTrials & uncuedTrials, 1, 'FluorDataField', 'ZS', 'sortValues', TE.lickLatency_us, 'zeroTimes', TE.Us, 'window', [-3 3]);
% subplot(1,6,6); phRasterFromTE(TE, rewardTrials & uncuedTrials, 2, 'FluorDataField', 'ZS', 'sortValues', TE.lickLatency_us, 'zeroTimes', TE.Us, 'window', [-3 3]);
% 
% 
% ensureFigure('Reward_vs_LickRate', 1);
% 
% % lickrates = sort(TE.licks_us.rate(rewardTrials & ~uncuedTrials));
% % subplot(1,6,1); plot(lickrates, 1:sum(rewardTrials & ~uncuedTrials), '-k'); set(gca, 'YDir', 'reverse');
% subplot(1,6,1); eventRasterFromTE(TE, rewardTrials & ~uncuedTrials, 'Port1In', 'trialNumbering', 'consecutive',...
%     'zeroTimes', TE.Outcome, 'window', xwindow, 'sortValues', TE.licks_us.rate);
% set(gca, 'XLim', xwindow);
% subplot(1,6,2); phRasterFromTE(TE, rewardTrials & ~uncuedTrials, 1, 'FluorDataField', 'ZS', 'sortValues', TE.licks_us.rate	, 'zeroTimes', TE.Us, 'window', [-3 3]);
% subplot(1,6,3); phRasterFromTE(TE, rewardTrials & ~uncuedTrials, 2, 'FluorDataField', 'ZS', 'sortValues', TE.licks_us.rate, 'zeroTimes', TE.Us, 'window', [-3 3]);
% lickrates = sort(TE.licks_us.rate(rewardTrials & uncuedTrials));
% % subplot(1,6,4); plot(lickrates, 1:sum(rewardTrials & uncuedTrials), '-k'); set(gca, 'YDir', 'reverse');
% % set(gca, 'XLim', xwindow);
% subplot(1,6,4); eventRasterFromTE(TE, rewardTrials & uncuedTrials, 'Port1In', 'trialNumbering', 'consecutive',...
%     'zeroTimes', TE.Outcome, 'window', xwindow, 'sortValues', TE.licks_us.rate);
% set(gca, 'XLim', xwindow);
% subplot(1,6,5); phRasterFromTE(TE, rewardTrials & uncuedTrials, 1, 'FluorDataField', 'ZS', 'sortValues', TE.licks_us.rate, 'zeroTimes', TE.Us, 'window', [-3 3]);
% subplot(1,6,6); phRasterFromTE(TE, rewardTrials & uncuedTrials, 2, 'FluorDataField', 'ZS', 'sortValues', TE.licks_us.rate, 'zeroTimes', TE.Us, 'window', [-3 3]);