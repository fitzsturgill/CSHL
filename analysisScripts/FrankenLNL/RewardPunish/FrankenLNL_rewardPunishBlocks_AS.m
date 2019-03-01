%% 1/21/19 FS

saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_FrankenLNL_4odors(sessions);

%%
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRCaMP(ch2)
channels=[]; dFFMode = {}; BL = {}; 
if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'expFit';
%     dFFMode{end+1} = 'simple';    
    BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'expFit';
    BL{end + 1} = [1 4];    
end

%% process photometry
% baseline by trial
 TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue2', 'channels', channels, 'baseline', BL);
% baseline expfit
TE.PhotometryExpFit = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'expFit', 'zeroField', 'Cue2', 'channels', channels, 'baseline', BL);
% TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'expFitBegin', 0.1,...
%     'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL, 'downsample', 305);
duration = length(TE.Photometry.xData) / TE.Photometry.sampleRate;




%% add wheel
duration = length(TE.Photometry.xData) / TE.Photometry.sampleRate;
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', duration, 'Fs', 20, 'startField', 'PreCsRecording', 'zeroField', 'Cue2');

%% add eye avg if desired/present
duration = length(TE.Photometry.xData) / TE.Photometry.sampleRate;
TE.eyeAvg = addCSVToTE(TE, 'duration', duration, 'zeroField', 'Outcome', 'folderPrefix', 'EyeAvg_', 'filePrefix', 'EyeAvg_');%, 'normMode', 'byTrial');

%% add pupil if desired/present
duration = length(TE.Photometry.xData) / TE.Photometry.sampleRate;
TE = addPupilometryToTE(TE, 'duration', duration, 'zeroField', 'Outcome', 'frameRate', 60, 'frameRateNew', 20, 'normMode', 'byTrial');
%%
% basepath = uigetdir;
basepath = 'Z:\SummaryAnalyses\Franken_LNL_aversive';
sep = strfind(TE.filename{1}, '_');
if length(unique(TE.filename)) > 1
    subjectName = TE.filename{1}(1:sep(2)-1);
else
    subjectName = TE.filename{1}(1:end-4);
end
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

%% truncate sessions
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'trialType', 1));

%% also reject shocks under 200uAmps if present
if isfield(TE, 'ShockCurrent')
    TE.reject = TE.reject | ((abs(TE.ShockCurrent) < 200) & (strcmp(TE.ReinforcementOutcome, 'Shock')));
end
%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%% cross sessions bleaching curve and exponential fits
for channel = channels
    figname = ['sessionBleach_Correction_ch' num2str(channel)];
    ensureFigure(figname, 1);
    plot(TE.PhotometryExpFit.data(channel).blF_raw, 'k'); hold on;
    plot(TE.PhotometryExpFit.data(channel).blF, 'r');
    if saveOn
        saveas(gcf, fullfile(savepath, [figname '.fig']));
        saveas(gcf, fullfile(savepath, [figname '.jpg']));
    end
    % cross trial bleaching fits for each session plotted as axis array
    if 1 %channel == 1
        figname = ['trialBleach_Correction_ch' num2str(channel)];
        f1 = ensureFigure(figname, 1);
        nSessions = size(TE.Photometry.bleachFit, 1);
        subA = ceil(sqrt(nSessions));
        for counter = 1:nSessions
            subplot(subA, subA, counter);
            plot(TE.Photometry.bleachFit(counter, channel).trialTemplateFullX, TE.Photometry.bleachFit(counter, channel).trialTemplateFull, 'g'); hold on;            
            plot(TE.Photometry.bleachFit(counter, channel).trialTemplateX, TE.Photometry.bleachFit(counter, channel).trialTemplate, 'b'); 
            plot(TE.Photometry.bleachFit(counter, channel).fitX, TE.Photometry.bleachFit(counter, channel).trialFit, 'r');
        %     title(num2str(counter));    
        end
        % average of all trials for this channel to eyeball correction
        figname2 = ['corrected_allTrials_ch' num2str(channel)]; 
        ensureFigure(figname2, 1);
        [ha, hl] = phPlotAverageFromTE(TE, 1:length(TE.filename), channel,...
    'FluorDataField', 'ZS', 'window', [0.1, max(TE.Photometry.xData) - min(TE.Photometry.xData)], 'zeroTimes', TE.Photometry.startTime); %high value, reward
        if saveOn
            saveas(f1, fullfile(savepath, figname), 'fig');
            saveas(f1, fullfile(savepath, figname), 'jpeg');            
            saveas(gcf, fullfile(savepath, [figname2 '.fig']));
            saveas(gcf, fullfile(savepath, [figname2 '.jpg']));            
        end
    end
end

%% 
frankenLNL_conditions;

%% raster plots
PhotometryField = 'Photometry';
trialNumbering = 'consecutive';
CLimFactor = 3;
window = [-4 6];
for channel = channels
    saveName = ['Reward_phRasters_ch' num2str(channel)];
    ensureFigure(saveName, 1);
    
    subplot(1,5,1);
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
        'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('cued reward');         
    addStimulusPatch(gca, [0 1]);
    subplot(1,5,2);
    phRasterFromTE(TE, trialsByType{1}, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('cued reward'); set(gca, 'XLim', window);
    subplot(1,5,3);
    phRasterFromTE(TE, trialsByType{2}, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('cued omit');    set(gca, 'XLim', window);
    subplot(1,5,4);
    phRasterFromTE(TE, trialsByType{5}, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('uncued reward');       set(gca, 'XLim', window);
    
    subplot(1,5,5);
    eventRasterFromTE(TE, trialsByType{5}, 'Port1In', 'trialNumbering', trialNumbering,...
        'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         title('uncued reward');     
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end
    
    saveName = ['Punish_phRasters_ch' num2str(channel)];
    ensureFigure(saveName, 1);
    
    subplot(1,3,1);
    phRasterFromTE(TE, trialsByType{3}, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window);
    xlabel('time from cue (s)');  title('cued punish');
    subplot(1,3,2);
    phRasterFromTE(TE, trialsByType{4}, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window);
    xlabel('time from cue (s)');  title('cued ommission');    
    subplot(1,3,3);
    phRasterFromTE(TE, trialsByType{6}, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window);
    xlabel('time from cue (s)');  title('uncued punish');       
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end    
end

%% blink rasters
saveName = 'Aversive_blinkRasters';
ensureFigure(saveName, 1); colormap jet;
% clims = [0 1.5];
cLimFactor = 4;
blinkData = TE.eyeAvg.EyeAvgNorm;
blinkDataMean = mean2(blinkData);
clims = [nanmean(blinkData(:)) - nanstd(blinkData(:)) * cLimFactor, nanmean(blinkData(:)) + nanstd(blinkData(:)) * cLimFactor];
xData = [TE.eyeAvg.xData(1) TE.eyeAvg.xData(end)];
subplot(1,3,1); imagesc('XData', xData, 'CData', blinkData(trialsByType{3}, :), clims); set(gca, 'XLim', xData);
title('cued punish');
subplot(1,3,2); imagesc('XData', xData, 'CData', blinkData(trialsByType{4}, :), clims); set(gca, 'XLim', xData);
title('cued ommission');
subplot(1,3,3); imagesc('XData', xData, 'CData', blinkData(trialsByType{6}, :), clims); set(gca, 'XLim', xData);
title('uncued punish');

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
[data, XData] = alignedDataWindow(TE, TE.pupil.pupDiameterNorm, true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window);
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
[data, XData] = alignedDataWindow(TE, TE.pupil.pupDiameterNorm, true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window);
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
[data, XData] = alignedDataWindow(TE, TE.pupil.(dataField), true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window);
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
[data, XData] = alignedDataWindow(TE.Wheel.data.V, true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window, 'Fs', 20);
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
[data, XData] = alignedDataWindow(TE, TE.pupil.pupDiameterNorm, true(nTrials, 1), 'startTimes', TE.Photometry.startTime, 'zeroTimes', TE.Cue2, 'window', window);
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

% Define the positions of different axes matrices on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    params = struct();
    matpos_title = [0 0 1 .1];
    matpos_avgs = [0 .1 1 0.85];    
    params.cellmargin = [.05 .05 0.05 0.05];    
    params.figmargin = [0.05 0.05 0.025 0.025];


    
    
    % -3 6
window = [-1 5];
PhotometryField = 'Photometry';
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);

% rewarding subset
showTheseR = find(Odor2Valve1Trials & rewardTrials);
[~, rixR] = sort(rand(size(showTheseR)));
% aversive subset
showTheseA = find(Odor2Valve2Trials & punishTrials);
[~, rixA] = sort(rand(size(showTheseA)));


for channel = channels
    saveName = ['Averages_ch' num2str(channel)];
    fig = ensureFigure(saveName, 1);
    mcPortraitFigSetup(fig);

    params.matpos = matpos_title;
    [ax, ~] = textAxes(fig, sprintf('%s , Channel %d', subjectName, channel), 14);
    setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);    
    
    try
        ax = subplot(3,2,1);
        plot(TE.Photometry.xData, TE.Photometry.data(channel).ZS(showTheseR(rixR(1:10)), :)', 'LineWidth', 0.25); set(gca, 'XLim', window);
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('Fluor. (ZS)'); title('rewarding');
        disp('hello');
    catch
    end
    
    try
        ax(2) = subplot(3,2,2);
        plot(TE.Photometry.xData, TE.Photometry.data(channel).ZS(showTheseA(rixA(1:10)), :)', 'LineWidth', 0.25); set(gca, 'XLim', window);
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('Fluor. (ZS)'); title('aversive');
    catch
    end    

    ax(3) = subplot(3,2,3);
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5]), channel, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k', 'c'}, 'PhotometryField', PhotometryField); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
    set(gca, 'XLim', window); ylabel('Fluor ZS');
    
    try
        ax(4) = subplot(3,2,4);
        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([3 4 6]), channel, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k', 'm'}, 'PhotometryField', PhotometryField); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        set(gca, 'XLim', window); ylabel('Fluor ZS');
    catch
    end
    
    try
        ax(5) = subplot(3,2,5);
        varargin = {'window', window, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k','c'}};
        axh = [];
        [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 2 5]), 'Port1In', varargin{:});  
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('licks (1/s)'); xlabel('time from cue(s)');  set(gca, 'XLim', window);
        legend(hl, {'cued reward', 'cued ommission', 'uncued reward'}, 'Box', 'off', 'Location', 'best'); 
    catch
    end
    
    try
        ax(6) = subplot(3,2,6);
        varargin = {'window', window, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k','m'}};
        axh = [];
        [ha, hl] = plotEventAverageFromTE(TE, trialsByType([3 4 6]), 'Port1In', varargin{:});  
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('licks (1/s)'); xlabel('time from cue(s)');  set(gca, 'XLim', window); legend(hl, {'cued punish', 'cued ommission', 'uncued punish'}, 'Box', 'off', 'Location', 'best');      
    catch
    end
    

    
    params.matpos = matpos_avgs;
    setaxesOnaxesmatrix(ax, 3, 2, 1:6, params, fig);     
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
[~, rix] = sort(rand(size(showThese)));
plot(TE.Photometry.xData, TE.Photometry.data(1).ZS(showThese(rix(1:20)), :)', 'k', 'LineWidth', 0.25); set(gca, 'XLim', [-2 4]);
% plot(TE.Photometry.xData, TE.Photometry.data(1).ZS(showThese(rix(1:10)), :)', 'LineWidth', 0.25); set(gca, 'XLim', [-2 5]);
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);
addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
xlabel('Time from odor (s)'); ylabel('Fluor. (ZS)'); title('ACh. sensor in BLA');
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end    
    
    
%% synchrony between left and right BLA
saveName = 'examples_synchrony_leftRightBLA';
ensureFigure(saveName, 1); 
showThese = find(Odor2Valve1Trials & rewardTrials);
[~, rix] = sort(rand(size(showThese)));
whichOne = 3;
hl = zeros(2,1);
hl(1) = plot(TE.Photometry.xData, TE.Photometry.data(1).ZS(showThese(rix(1)), :)', 'g', 'LineWidth', 1); set(gca, 'XLim', [-2 4]); hold on;
hl(2) = plot(TE.Photometry.xData, TE.Photometry.data(2).ZS(showThese(rix(1)), :)', 'r', 'LineWidth', 1); set(gca, 'XLim', [-2 4]);
addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [1.9 2.1]);
legend(hl, {'Flat', 'Tapered'}, 'Box', 'off', 'Location', 'northwest');
xlabel('time from cue (s)'); ylabel('Fluor. (ZS)');


if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end    

    
