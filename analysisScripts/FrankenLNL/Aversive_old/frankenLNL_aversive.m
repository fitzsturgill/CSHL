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
    dFFMode{end+1} = 'simple';
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

%% add pupilometry
% duration = 8.3;
TE = addPupilometryToTE(TE, 'duration', duration, 'zeroField', 'Cue2',  'frameRate', 60, 'frameRateNew', 20, 'normMode', 'byTrial');


%% add whisking
TE.Whisk = addWhiskingToTE(TE, 'duration', duration);

%% add wheel
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', duration, 'Fs', 20, 'startField', 'Start');

    sep = strfind(TE.filename{1}, '.');
    subjectName = TE.filename{1}(1:sep(1)-1);
%%
% basepath = uigetdir;
basepath = 'Z:\SummaryAnalyses\Franken_LNL_aversive';
sep = strfind(TE.filename{1}, '_');
if length(sessions) > 1
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
%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%% cross sessions bleaching curve and exponential fits
for channel = channels
    figname = ['sessionBleach_Correction_ch' num2str(channel)];
    ensureFigure(figname, 1);
    plot(TE.Photometry.data(channel).blF_raw, 'k'); hold on;
    plot(TE.Photometry.data(channel).blF, 'r');
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
CLimFactor = 3;
for channel = channels
    saveName = ['Aversive_phRasters_ch' num2str(channel)];
    ensureFigure(saveName, 1);
    
    subplot(1,3,1);
    phRasterFromTE(TE, ~uncuedTrials & punishTrials, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-4 6]);
    xlabel('time from cue (s)');  title('cued punish');
    subplot(1,3,2);
    phRasterFromTE(TE, ~uncuedTrials & neutralTrials, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-4 6]);
    xlabel('time from cue (s)');  title('cued ommission');    
    subplot(1,3,3);
    phRasterFromTE(TE, uncuedTrials & punishTrials, channel, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-4 6]);
    xlabel('time from cue (s)');  title('uncued punish');       
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end
end

%% whisk rasters
saveName = 'Aversive_whiskRasters';
ensureFigure(saveName, 1);
cLimFactor = 2;
clim = [mean2(TE.Whisk.whiskNorm) - cLimFactor * std2(TE.Whisk.whiskNorm) mean2(TE.Whisk.whiskNorm) + cLimFactor * std2(TE.Whisk.whiskNorm)];
subplot(1,3,1); imagesc(TE.Whisk.whiskNorm(~uncuedTrials & punishTrials, :), clim);
title('cued punish');
subplot(1,3,2); imagesc(TE.Whisk.whiskNorm(~uncuedTrials & neutralTrials, :), clim);
title('cued ommission');
subplot(1,3,3); imagesc(TE.Whisk.whiskNorm(uncuedPunish, :), clim);
title('uncued punish');

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end

%% wheel rasters and averages
saveName = 'Aversive_wheelRasters';
ensureFigure(saveName, 1);
subplot(2,3,1); imagesc(TE.Wheel.data.V(~uncuedTrials & punishTrials, :));
title('cued punish');
subplot(2,3,2); imagesc(TE.Wheel.data.V(~uncuedTrials & neutralTrials, :));
title('cued ommission');
subplot(2,3,3); imagesc(TE.Wheel.data.V(uncuedPunish, :));
title('uncued punish');
subplot(2,3,4); plot(nanmean(TE.Wheel.data.V(~uncuedTrials & punishTrials, :)));
subplot(2,3,5); plot(nanmean(TE.Wheel.data.V(~uncuedTrials & neutralTrials, :)));
subplot(2,3,6); plot(nanmean(TE.Wheel.data.V(uncuedPunish, :)));

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end


%% averages
saveName = 'Aversive_phAverages';
ensureFigure(saveName, 1);
trialSets = {uncuedTrials & punishTrials, ~uncuedTrials & neutralTrials, ~uncuedTrials & punishTrials};
titles = {'Ach.', 'Dop.'};
legend_text = {'uncued punish', 'cued ommission', 'cued punish'};
for channel = channels
    subplot(1,2,channel);
    [ha, hl] = phPlotAverageFromTE(TE, trialSets, channel, 'FluorDataField', 'ZS', 'window', [-4 6], 'linespec', {'m', 'k', 'r'}); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    xlabel('time from cue (s)');set(gca, 'XLim', [-4 6]); title(titles{channel}); legend(hl, legend_text, 'Box', 'off', 'Location', 'best');
end

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end
    
%% blink rasters
saveName = 'Aversive_blinkRasters';
ensureFigure(saveName, 1); colormap jet;
% clims = [0 1.5];
cLimFactor = 2;
blinkData = TE.pupil.eyeAreaNorm;
blinkDataMean = mean2(blinkData);
clims = [nanmean(blinkData(:)) - nanstd(blinkData(:)) * cLimFactor, nanmean(blinkData(:)) + nanstd(blinkData(:)) * cLimFactor];

subplot(1,3,1); imagesc('XData', [ -4.000 7], 'CData', blinkData(uncuedTrials & punishTrials, :), clims);
title('uncued punish');
subplot(1,3,2); imagesc('XData', [ -4 7], 'CData', blinkData(~uncuedTrials & neutralTrials, :), clims);
title('cued ommission');
subplot(1,3,3); imagesc('XData', [-4 7], 'CData', blinkData(~uncuedTrials & punishTrials, :), clims);
title('cued punish');

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end
%%
ensureFigure('test', 1);
clims = [0 1.5];
blinkData = TE.pupil.eyeAreaNorm;
% blinkData(isnan(blinkData)) = 2;
subplot(1,2,1); imagesc('XData', [ -1 1], 'CData', blinkData(:, bpX2pnt(-1, 20, -4):bpX2pnt(1, 20, -4)), clims);
subplot(1,2,2); imagesc('XData', [ -4 7], 'CData', blinkData(:, :), clims);


%% blink averages
saveName = 'Aversive_blinkAverages';
ensureFigure(saveName, 1);
xData = TE.pupil.xData;
blinkData = TE.pupil.eyeAreaNorm;
blinkData(isnan(blinkData)) = 0;
axes; hold on; grid on;
plot(xData, mean(blinkData(uncuedTrials & punishTrials, :)), 'r');
plot(xData, mean(blinkData(~uncuedTrials & neutralTrials, :)), 'k');
plot(xData, mean(blinkData(~uncuedTrials & punishTrials, :)), 'm-*');
title('frame avg (goes up as eye closes)');
xlabel('time from cue');
set(gca, 'XLim', [-1 2]);
legend('uncued', 'omit', 'cued'); 

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end
