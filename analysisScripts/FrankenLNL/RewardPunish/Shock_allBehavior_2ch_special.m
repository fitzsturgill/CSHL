

%% cued Reward, cued Shock, all behavior
PhotometryField = 'Photometry';
trialNumbering = 'global';
CLimFactor = 3;
window = [-4 6];

saveName = 'cued_reinforcment_allBehavior_shock_special';
ensureFigure(saveName, 1);

subplot(1,8,1);
eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
    'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('Cued Reward:');         
addStimulusPatch(gca, [0 2]);

subplot(1,8,2);
phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('HDB/ChAT'); set(gca, 'XLim', window, 'YTick', []);

subplot(1,8,3);
phRasterFromTE(TE, trialsByType{1}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('VTA/DAT'); set(gca, 'XLim', window, 'YTick', []);



subplot(1,8,4);
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('Cued Shock: HDB/ChAT'); set(gca, 'XLim', window, 'YTick', []);

subplot(1,8,5);
phRasterFromTE(TE, trialsByType{3}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('VTA/DAT'); set(gca, 'XLim', window, 'YTick', []);

subplot(1,8,6);
alignedDataRaster(TE.pupil.pupDiameterNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Pupil diameter'); set(gca, 'XLim', window, 'YTick', []);

subplot(1,8,7);
alignedDataRaster(TE.pupil.frameAvgNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
title('Eye closure'); set(gca, 'XLim', window, 'YTick', []);

subplot(1,8,8);
alignedDataRaster(TE.Wheel.data.V, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.Wheel.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.Wheel.Fs, 'CLim', [0 20]);
title('Velocity'); set(gca, 'XLim', window, 'YTick', []);

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end




%% pupil rasters
window = [-4 6];
saveName = 'Aversive_pupilRasters';
ensureFigure(saveName, 1); colormap jet;
% first just get all the data
nTrials = length(TE.filename);
[data, XData] = alignedDataWindow(TE.pupil.pupDiameterNorm, true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window);
% determine the clims
cmean = nanmean(nanmean(data(:, 1:20*4), 2), 1);
cstd = nanmean(nanstd(data(:, 1:20*4), 0, 2), 1);
cLimFactor = 4;
clims = [cmean - cLimFactor * cstd cmean + cLimFactor * cstd];
dummy = NaN(size(data));
subplot(1,3,1);
thisData = dummy;
thisData(trialsByType{3}, :) = data(trialsByType{3}, :);
imagesc(window, [1 nTrials], thisData, clims);
title('cued shock'); ylabel('trial #');
subplot(1,3,2);
thisData = dummy;
thisData(trialsByType{4}, :) = data(trialsByType{4}, :);
imagesc(window, [1 nTrials], thisData, clims);
title('omission'); xlabel('time from odor (s)');
subplot(1,3,3);
thisData = dummy;
thisData(trialsByType{6}, :) = data(trialsByType{6}, :);
imagesc(window, [1 nTrials], thisData, clims);
title('uncued shock');

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end


saveName = 'Appetitive_pupilRasters';
ensureFigure(saveName, 1); colormap jet;
% first just get all the data
nTrials = length(TE.filename);
[data, XData] = alignedDataWindow( TE.pupil.pupDiameterNorm, true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window);
% determine the clims
cmean = nanmean(nanmean(data(:, 1:20*4), 2), 1);
cstd = nanmean(nanstd(data(:, 1:20*4), 0, 2), 1);
cLimFactor = 4;
clims = [cmean - cLimFactor * cstd cmean + cLimFactor * cstd];
dummy = NaN(size(data));
subplot(1,3,1);
thisData = dummy;
thisData(trialsByType{1}, :) = data(trialsByType{1}, :);
imagesc(window, [1 nTrials], thisData, clims);
title('cued reward'); ylabel('trial #');
subplot(1,3,2);
thisData = dummy;
thisData(trialsByType{2}, :) = data(trialsByType{2}, :);
imagesc(window, [1 nTrials], thisData, clims);
title('omission'); xlabel('time from odor (s)');
subplot(1,3,3);
thisData = dummy;
thisData(trialsByType{5}, :) = data(trialsByType{5}, :);
imagesc(window, [1 nTrials], thisData, clims);
title('uncued reward');

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end

%% pupil rasters, consecutive numbering
window = [-4 6];
saveName = 'Aversive_pupilRasters';
ensureFigure(saveName, 1); colormap parula;
dataField = 'pupDiameterNorm';
% dataField = 'frameAvgNorm';
% first just get all the data
nTrials = length(TE.filename);
[data, XData] = alignedDataWindow( TE.pupil.(dataField), true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window);
% determine the clims
cmean = nanmean(nanmean(data(:, 1:20*4), 2), 1);
cstd = nanmean(nanstd(data(:, 1:20*4), 0, 2), 1);
cLimFactor = 4;
clims = [cmean - cLimFactor * cstd cmean + cLimFactor * cstd];

subplot(1,3,1);

thisData = data(trialsByType{3}, :);
imagesc(window, [1 sum(trialsByType{3})], thisData, clims);
sessionBreaks = find(diff(TE.sessionIndex(trialsByType{3})))';            
line(repmat(window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
title('cued shock'); ylabel('trial #');
subplot(1,3,2);
thisData = data(trialsByType{4}, :);
imagesc(window, [1 sum(trialsByType{4})], thisData, clims);
sessionBreaks = find(diff(TE.sessionIndex(trialsByType{4})))';            
line(repmat(window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
title('omission'); xlabel('time from odor (s)');
subplot(1,3,3);

thisData = data(trialsByType{6}, :);
imagesc(window, [1 sum(trialsByType{6})], thisData, clims);
sessionBreaks = find(diff(TE.sessionIndex(trialsByType{6})))';            
line(repmat(window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
title('uncued shock');

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end


saveName = 'Appetitive_pupilRasters';
ensureFigure(saveName, 1); colormap parula;
% first just get all the data


subplot(1,3,1);
thisData = data(trialsByType{1}, :);
imagesc(window, [1 nTrials], thisData, clims);
title('cued reward'); ylabel('trial #');
subplot(1,3,2);
thisData = data(trialsByType{2}, :);
imagesc(window, [1 nTrials], thisData, clims);
title('omission'); xlabel('time from odor (s)');
subplot(1,3,3);
thisData = data(trialsByType{5}, :);
imagesc(window, [1 nTrials], thisData, clims);
title('uncued reward');

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end

%% aversive wheel raster

window = [-4 6];
saveName = 'Aversive_wheelRaster';
ensureFigure(saveName, 1); colormap parula;
% first just get all the data
nTrials = length(TE.filename);
[data, XData] = alignedDataWindow( TE.Wheel.data.V, true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window);
% determine the clims
cmean = nanmean(nanmean(data(:, 1:20*4), 2), 1);
cstd = nanmean(nanstd(data(:, 1:20*4), 0, 2), 1);
cLimFactor = 20;
clims = [cmean - cLimFactor * cstd cmean + cLimFactor * cstd];
clims = [0 20];

subplot(1,3,1);
thisData = data(trialsByType{3}, :);
imagesc(window, [1 sum(trialsByType{3})], thisData, clims);
sessionBreaks = find(diff(TE.sessionIndex(trialsByType{3})))';            
line(repmat(window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
title('cued shock'); ylabel('trial #');

subplot(1,3,2);
thisData = data(trialsByType{4}, :);
imagesc(window, [1 sum(trialsByType{4})], thisData, clims);
sessionBreaks = find(diff(TE.sessionIndex(trialsByType{4})))';            
line(repmat(window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
title('omission'); xlabel('time from odor (s)');

subplot(1,3,3);
thisData = data(trialsByType{6}, :);
imagesc(window, [1 sum(trialsByType{6})], thisData, clims);
sessionBreaks = find(diff(TE.sessionIndex(trialsByType{6})))';            
line(repmat(window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
title('uncued shock');

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end

%% pupil rasters, consecutive numbering
window = [-4 6];
saveName = 'Aversive_pupilRasters';
ensureFigure(saveName, 1); colormap parula;
% first just get all the data
nTrials = length(TE.filename);
[data, XData] = alignedDataWindow( TE.pupil.pupDiameterNorm, true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window);
% determine the clims
cmean = nanmean(nanmean(data(:, 1:20*4), 2), 1);
cstd = nanmean(nanstd(data(:, 1:20*4), 0, 2), 1);
cLimFactor = 4;
clims = [cmean - cLimFactor * cstd cmean + cLimFactor * cstd];

subplot(1,3,1);
thisData = data(trialsByType{3}, :);
imagesc(window, [1 sum(trialsByType{3})], thisData, clims);
sessionBreaks = find(diff(TE.sessionIndex(trialsByType{3})))';            
line(repmat(window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
title('cued shock'); ylabel('trial #');

subplot(1,3,2);
thisData = data(trialsByType{4}, :);
imagesc(window, [1 sum(trialsByType{4})], thisData, clims);
sessionBreaks = find(diff(TE.sessionIndex(trialsByType{4})))';            
line(repmat(window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
title('omission'); xlabel('time from odor (s)');

subplot(1,3,3);
thisData = data(trialsByType{6}, :);
imagesc(window, [1 sum(trialsByType{6})], thisData, clims);
sessionBreaks = find(diff(TE.sessionIndex(trialsByType{6})))';            
line(repmat(window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
title('uncued shock');

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end

%% averages
% -3 6
window = [-1 4];
for channel = channels
    saveName = ['Averages_ch' num2str(channel)];
    ensureFigure(saveName, 1);
    
    subplot(2,2,1);
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5]), channel, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k', 'c'}); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    set(gca, 'XLim', window); title('rewarding'); ylabel('Fluor ZS');
    
    try
        subplot(2,2,2);
        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([3 4 6]), channel, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k', 'm'}); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        set(gca, 'XLim', window); title('aversive'); ylabel('Fluor ZS');
    catch
    end
    
    try
        subplot(2,2,3);
        varargin = {'window', window, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k','c'}};
        axh = [];
        [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 2 5]), 'Port1In', varargin{:});  
        ylabel('licks (1/s)'); xlabel('time from cue(s)');  set(gca, 'XLim', window);
        legend(hl, {'cued reward', 'cued ommission', 'uncued reward'}, 'Box', 'off', 'Location', 'best'); 
    catch
    end
    
    try
        subplot(2,2,4);
        varargin = {'window', window, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k','m'}};
        axh = [];
        [ha, hl] = plotEventAverageFromTE(TE, trialsByType([3 4 6]), 'Port1In', varargin{:});  
        ylabel('licks (1/s)'); xlabel('time from cue(s)');  set(gca, 'XLim', window); legend(hl, {'cued punish', 'cued ommission', 'uncued punish'}, 'Box', 'off', 'Location', 'best');      
    catch
    end
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end    
end

%% blink averages
saveName = 'Aversive_blinkAverages';
ensureFigure(saveName, 1);
xData = TE.eyeAvg.xData;
blinkData = TE.eyeAvg.EyeAvgNorm;
blinkData(isnan(blinkData)) = 0;
axes; hold on; grid on;
boundedline(xData(:), [mean(blinkData(trialsByType{3}, :)); mean(blinkData(trialsByType{4}, :)); mean(blinkData(trialsByType{6}, :))]',...
    permute([std(blinkData(trialsByType{3}, :)) ./ sqrt(sum(trialsByType{3})); std(blinkData(trialsByType{4}, :)) ./ sqrt(sum(trialsByType{3})); std(blinkData(trialsByType{6}, :)) ./ sqrt(sum(trialsByType{3}))], [2 3 1]),...
    'cmap', [1 0 0; 0 0 0; 1 0 1]);
% plot(xData, mean(blinkData(trialsByType{3}, :)), 'r');
% plot(xData, mean(blinkData(trialsByType{4}, :)), 'k');
% plot(xData, mean(blinkData(trialsByType{6}, :)), 'm-*');
title('frame avg (goes up as eye closes)');
xlabel('time from punishment');
set(gca, 'XLim', [-4 3]);
legend('cued', 'omit', 'uncued'); 

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end

%% wheel averages
saveName = 'Aversive_wheelAverages';
ensureFigure(saveName, 1);
xData = TE.Wheel.xData;
wheelData = TE.Wheel.data.V;
wheelData(isnan(wheelData)) = 0;
axes; hold on; grid on;
boundedline(xData(:), [nanmean(wheelData(trialsByType{3}, :)); nanmean(wheelData(trialsByType{4}, :)); nanmean(wheelData(trialsByType{6}, :))]',...
    permute([nanstd(wheelData(trialsByType{3}, :)) ./ sqrt(sum(isfinite(wheelData(trialsByType{3},:)), 1));...
    nanstd(wheelData(trialsByType{4}, :)) ./ sqrt(sum(isfinite(wheelData(trialsByType{4},:)), 1));...
    nanstd(wheelData(trialsByType{6}, :)) ./ sqrt(sum(isfinite(wheelData(trialsByType{6},:)), 1))], [2 3 1]),...
    'cmap', [1 0 0; 0 0 0; 1 0 1]);

title('wheel avg');
xlabel('time from cue'); ylabel('Velocity');
set(gca, 'XLim', [-4 7]);
legend('cued', 'omit', 'uncued'); 

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end

%% pupil averages
saveName = 'Aversive_pupilAverages';
ensureFigure(saveName, 1);
xData = TE.pupil.xData;
pupData = TE.pupil.pupDiameterNorm;
axes; hold on; grid on;
% boundedline(xData(:), [mean(pupData(trialsByType{3}, :)); mean(pupData(trialsByType{4}, :)); mean(pupData(trialsByType{6}, :))]',...
%     permute([std(pupData(trialsByType{3}, :)) ./ sqrt(sum(trialsByType{3})); std(pupData(trialsByType{4}, :)) ./ sqrt(sum(trialsByType{4})); std(pupData(trialsByType{6}, :)) ./ sqrt(sum(trialsByType{6}))], [2 3 1]),...
%     'cmap', [1 0 0; 0 0 0; 1 0 1]);

boundedline(xData(:), [nanmean(pupData(trialsByType{3}, :)); nanmean(pupData(trialsByType{4}, :)); nanmean(pupData(trialsByType{6}, :))]',...
    permute([nanstd(pupData(trialsByType{3}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{3},:)), 1));...
    nanstd(pupData(trialsByType{4}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{4},:)), 1));...
    nanstd(pupData(trialsByType{6}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{6},:)), 1))], [2 3 1]),...
    'cmap', [1 0 0; 0 0 0; 1 0 1]);

title('pupil avg');
xlabel('time from punishment'); ylabel('Pup diameter norm');
set(gca, 'XLim', [-4 3]);
legend('cued', 'omit', 'uncued'); 

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end

saveName = 'Appetitive_pupilAverages';
ensureFigure(saveName, 1);
xData = TE.pupil.xData;
pupData = TE.pupil.pupDiameterNorm;
axes; hold on; grid on;

boundedline(xData(:), [nanmean(pupData(trialsByType{1}, :)); nanmean(pupData(trialsByType{2}, :)); nanmean(pupData(trialsByType{5}, :))]',...
    permute([nanstd(pupData(trialsByType{1}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{1},:)), 1));...
    nanstd(pupData(trialsByType{2}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{2},:)), 1));...
    nanstd(pupData(trialsByType{5}, :)) ./ sqrt(sum(isfinite(pupData(trialsByType{5},:)), 1))], [2 3 1]),...
    'cmap', [0 0 1; 0 0 0; 0 1 1]);

title('pupil avg');
xlabel('time from reward'); ylabel('Pup diameter norm');
set(gca, 'XLim', [-4 3]);
legend('cued', 'omit', 'uncued'); 

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end



%%
saveName = 'examples_cuedReward_ch1';
showThese = find(Odor2Valve1Trials & rewardTrials);
ensureFigure(saveName);
plot(TE.Photometry.xData, TE.Photometry.data(1).ZS(showThese(1:20), :)', 'k', 'LineWidth', 0.25); set(gca, 'XLim', [-2 4]);
figure; plot(TE.Photometry.xData, TE.Photometry.data(1).ZS(showThese(1:20), :)', 'k', 'LineWidth', 0.25); set(gca, 'XLim', [-2 4]);

addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [1.9 2.1]);
xlabel('Time from odor (s)'); ylabel('Fluor. (ZS)'); title('ACh. sensor in BLA');
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end    
    
    
%% synchrony between left and right BLA
saveName = 'examples_synchrony_leftRightBLA';
showThese = find(Odor2Valve1Trials & rewardTrials);
ensureFigure(saveName, 1); 
whichOne = 3;
hl = zeros(2,1);
hl(1) = plot(TE.Photometry.xData, TE.Photometry.data(1).ZS(showThese(whichOne), :)', 'g', 'LineWidth', 1); set(gca, 'XLim', [-2 4]); hold on;
hl(2) = plot(TE.Photometry.xData, TE.Photometry.data(2).ZS(showThese(whichOne), :)', 'r', 'LineWidth', 1); set(gca, 'XLim', [-2 4]);
addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [1.9 2.1]);
legend(hl, {'Flat', 'Tapered'}, 'Box', 'off', 'Location', 'northwest');
xlabel('time from cue (s)'); ylabel('Fluor. (ZS)');


if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end    

    
