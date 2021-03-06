% 4/10/17  Analysis script for pavlovian reversals using LickNoLick_Odor_V2
% protocol

saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_LNL_odor_V2(sessions);
%%
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRCaMP (ch2)
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

%% baseline by trial
 TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL);
%% baseline expfit
TE.PhotometryExpFit = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', {'expFit', 'expFit'}, 'zeroField', 'Cue', 'channels', channels, 'baseline', BL,...
    'tau', 2);
% TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'expFitBegin', 0.1,...
%     'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL, 'downsample', 305);
%%
% if you are reloading TE do this:
channels = [1 2];
%% extract peak trial dFF responses to cues and reinforcement and lick counts
% csWindow = cellfun(&#64;(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);

csWindow = zeros(nTrials, 2);
csWindow(:,2) = cellfun(@(x,y,z) max(x(end), y(end)) - z(1), TE.AnswerLick, TE.AnswerNoLick, TE.Cue); 
% max 1) to select AnswerNoLick time stampfor no lick trials (unused state contains NaN)
% 2) to select AnswerLick time stamp for lick trials (AnswerLick follows
% AnswerNoLick state)

% set time windows to compute cue/conditioned stimulus (Cs) response for
% each trial
ch1CsWindow = csWindow;
ch2CsWindow = csWindow;

usWindow = [0 0.75];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'

% percentile value for peak estimations
percentValue = 0.8;

% estimate respones for different events in each trial for photometry
% (bpCalcPeak_dFF) and for licking (countEventFromTE)
TE.csLicks = countEventFromTE(TE, 'Port1In', [0 3], TE.Cue);
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);



for channel = channels
    TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'mean', 'phField', 'ZS');
    if channel == 1
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
    elseif channel == 2
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');        
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');  
    end
end

% TE.wheel_cs = mean(TE.Wheel.data.V(:,bpX2pnt(ch1CsWindow(1) + pupLag, 20, -4):bpX2pnt(ch1CsWindow(2) + pupLag, 20, -4)), 2);


%% add pupilometry
TE = addPupilometryToTE(TE, 'duration', 11, 'zeroField', 'Cue',  'frameRate', 60, 'frameRateNew', 20);
pupLag = 0;
TE.pupil_cs = bpCalcPeak_Pupil(TE.pupil, 'zeroTimes', TE.Cue, 'window', csWindow);
TE.pupil_us = bpCalcPeak_Pupil(TE.pupil, 'zeroTimes', TE.Us, 'window', [0.5 1.5]);
TE.pupil_baseline = bpCalcPeak_Pupil(TE.pupil, 'zeroTimes', TE.PreCsRecording, 'window', [0 4]);

%% add whisking
TE.Whisk = addWhiskingToTE(TE);
TE.whisk_cs = bpCalcPeak_Whisk(TE.Whisk, 'zeroTimes', TE.Cue, 'window', csWindow);
TE.whisk_baseline = bpCalcPeak_Whisk(TE.Whisk, 'zeroTimes', TE.PreCsRecording, 'window', [0 4]);

%% add wheel
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', 11, 'Fs', 20, 'startField', 'Start');
TE.wheel_baseline = mean(TE.Wheel.data.V(:,bpX2pnt(-4, 20, -4):bpX2pnt(0, 20, -4)), 2);
TE.wheel_cs = mean(TE.Wheel.data.V(:,bpX2pnt(csWindow(1), 20, -4):bpX2pnt(csWindow(2), 20, -4)), 2);

%%
% basepath = uigetdir; % prompts windows/mac osx to give you a location to save
basepath = 'Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk\';
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

%%
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
%     if channel == 1
        figname = ['trialBleach_Correction_ch' num2str(channel)];
        ensureFigure(figname, 1);
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
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.fig'));
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.jpg'));
            saveas(gcf, fullfile(savepath, [figname2 '.fig']));
            saveas(gcf, fullfile(savepath, [figname2 '.jpg']));            
        end
%     end
end
%% generate trial lookups for different combinations of conditions
LNL_conditions;


%% photometry averages, zscored
%     ylim = [-2 8];
    fdField = 'ZS';
    saveName = sprintf('%s_phAvgs_%s', subjectName, fdField);  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);

    pm = [2 2];
    
    % - 6 0 4
    if ismember(1, channels)
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, neutralTrials, uncuedReward}, 1,...
            'FluorDataField', fdField, 'window', [1, 7], 'linespec', {'b','k','c'}); %high value, reward
        legend(hl, {'rew', 'neutral', 'uncued rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel(sprintf('BF %s', fdField)); textBox(subjectName);%set(gca, 'YLim', ylim);
    end
    
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, neutralTrials, uncuedReward}, 2,...
            'FluorDataField', fdField, 'window', [1, 7], 'linespec', {'b','k','c'}); %high value, reward
        legend(hl, {'rew', 'neutral', 'uncued rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        ylabel(sprintf('VTA %s', fdField)); xlabel('time from cue (s)'); %set(gca, 'YLim', ylim);
    end
    
    % - 6 0 4
    if ismember(1, channels)    
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 1,...
        'FluorDataField', fdField, 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('CS+, outcomes'); set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);
    end
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 2,...
            'FluorDataField', fdField, 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        xlabel('time from cue (s)');     set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);
    end
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end    
%% pupil averages per session, compare csPlus with uncued to check that trial/file alignment is correct given issues with bonsai
nSessions = max(TE.sessionIndex);
ensureFigure('pupil_check_sessions_avg', 1);
nax = ceil(sqrt(nSessions));
for counter = 1:nSessions
    subplot(nax,nax,counter); 
        plotPupilAverageFromTE(TE, {uncuedTrials & TE.sessionIndex == counter, csPlusTrials & TE.sessionIndex == counter});
        set(gca, 'XLim', [-2 6]);
end

%% 
saveName = ('Pupil_averages');
ensureFigure(saveName, 1); 
% condition on behavior and cue condition
subplot(2,2,1);[~, hl] = plotPupilAverageFromTE(TE, {uncuedTrials, csPlusTrials & hitTrials, csMinusTrials & CRTrials}, 'window', [-2 6]);
legend(hl, {'uncued', 'cs+, hit', 'cs-, CR'}, 'Box', 'off', 'Location', 'northwest');
title('Cue condition'); xlabel('time from cue (s)');
set(gca, 'XLim', [-2 6]);
% condition on behavior only for csPlus
subplot(2,2,2);[~, hl] = plotPupilAverageFromTE(TE, {hitTrials & csPlusTrials, missTrials & csPlusTrials}, 'window', [-2 6]);
set(gca, 'XLim', [-2 6]);
title('Cs+ by behavioral response'); xlabel('time from cue (s)');
legend(hl, {'cs+, hit', 'cs+, miss'}, 'Box', 'off', 'Location', 'northwest');
% reward, expected vs unexpected
subplot(2,2,3);[~, hl] = plotPupilAverageFromTE(TE, {csPlusTrials & hitTrials & rewardTrials, (csMinusTrials & CRTrials & rewardTrials) | uncuedReward}, 'window', [-2 6]);
set(gca, 'XLim', [-2 6]);
legend(hl, {'expected reward', 'unexpected reward'}, 'Box', 'off', 'Location', 'northwest');
title('Reward by cue and behavior'); xlabel('time from cue (s)');
formatFigure('scaleFactor', 4)
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
end
%% all behavior CsPlus
    fdField = 'diff';
    saveName = 'all_behavior_CsPlus';
    ensureFigure(saveName, 1);
    
    reversals = find(diff(TE.BlockNumber(csPlusTrials, :))) + 1;
    sessionChanges = find(diff(TE.sessionIndex(csPlusTrials, :))) + 1;

    
    
    subplot(1,6,1);    
    image(TE.Wheel.data.V(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 3, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 3]); 
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1);
    title('Velocity');
    ylabel('trial number');

    subplot(1,6,2);
    title('pupil');
    try
        image(TE.pupil.pupDiameterNorm(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
        set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
        line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
        colormap('parula');  
        title('Pupil Diameter');    
    catch
    end
    subplot(1,6,3);
    imagesc(TE.Whisk.whiskNorm(csPlusTrials, :), 'XData', [-4 7], [0 2])
    
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('whisking');
    
    subplot(1,6,4);
    eventRasterFromTE(TE, csPlusTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7]);
    title('licking');
    
    
    subplot(1,6,5); phRasterFromTE(TE, csPlusTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 1, 'FluorDataField', fdField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,6,6); phRasterFromTE(TE, csPlusTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 1, 'FluorDataField', fdField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('DAT');    
    axs = findobj(gcf, 'Type', 'axes');

    set(axs, 'FontSize', 16);
    set(axs(2:end), 'YTick', []);
    set(gcf, 'Position', [1 1 1920 1004]);
    saveas(gcf, fullfile(savepath, saveName), 'fig'); 
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');


%% all behavior CsMinus
    saveName = 'all_behavior_CsMinus';
    ensureFigure(saveName, 1);
    
    reversals = find(diff(TE.BlockNumber(csMinusTrials, :))) + 1;
    sessionChanges = find(diff(TE.sessionIndex(csMinusTrials, :))) + 1;
    
    subplot(1,6,1);    
    image(TE.Wheel.data.V(csMinusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines      
    title('Velocity');
    ylabel('trial number');

    subplot(1,6,2);
    title('pupil');
    try
        image(TE.pupil.pupDiameterNorm(csMinusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
        set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
        line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
        colormap('parula');  
        title('Pupil Diameter');    
    catch
    end
    subplot(1,6,3);
    imagesc(TE.Whisk.whiskNorm(csMinusTrials, :), 'XData', [-4 7], [0 2])
    
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('whisking');
    
    subplot(1,6,4);
    eventRasterFromTE(TE, csMinusTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7]);
    title('licking');
    
    
    subplot(1,6,5); phRasterFromTE(TE, csMinusTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    
    subplot(1,6,6); phRasterFromTE(TE, csMinusTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('DAT');    
    axs = findobj(gcf, 'Type', 'axes');

set(axs, 'FontSize', 16);
set(axs(2:end), 'YTick', []);
set(gcf, 'Position', [1 1 1920 1004]);
saveas(gcf, fullfile(savepath, 'allBehavior_whisk_csMinus'), 'fig'); 
saveas(gcf, fullfile(savepath, 'allBehavior_whisk_csMinus'), 'jpeg');



%% compile data into nReversals x nTrials arrays

    dataToPull = {...
        'csLicks', TE.csLicks.rate,...
        'usLicks', TE.usLicks.rate,...
        'phPeakMean_cs_ch1', TE.phPeakMean_cs(1).data,...
        'phBaseline_ch1', TE.phPeakMean_baseline(1).data,...
        'phPeakPercentile_cs_ch1', TE.phPeakPercentile_cs(1).data,...
        'phPeakMean_us_ch1', TE.phPeakMean_us(1).data,...
        'phPeakPercentile_us_ch1', TE.phPeakPercentile_us(1).data,...
        'trialOutcome', TE.trialOutcome,...
        'trialType', TE.trialType,...
        'trialNumber', TE.trialNumber,...
        'filename', TE.filename,...
        'ReinforcementOutcome', TE.ReinforcementOutcome,...
        'OdorValveIndex', TE.OdorValveIndex,...
        'csLicksROC', TE.AnswerLicksROC,...
        'whisk_cs', TE.whisk_cs.data,...
        'whisk_baseline', TE.whisk_baseline.data,...
        'pupil_cs', TE.pupil_cs.data,...
        'pupil_baseline', TE.pupil_baseline.data,...  
        'wheel_cs', TE.wheel_cs,...
        'wheel_baseline', TE.wheel_baseline,...
        };
    
    if ismember(2, channels)
        dataToPull = [dataToPull...
            {...
        'phPeakMean_cs_ch2', TE.phPeakMean_cs(2).data,...
        'phBaseline_ch2', TE.phPeakMean_baseline(2).data,...        
        'phPeakPercentile_cs_ch2', TE.phPeakPercentile_cs(2).data,...
        'phPeakMean_us_ch2', TE.phPeakMean_us(2).data,...
        'phPeakPercentile_us_ch2', TE.phPeakPercentile_us(2).data}...
        ];
    end
    RE = struct();
    RE.csPlus = extractReversalsFromTE(TE, csPlusTrials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    RE.csMinus = extractReversalsFromTE(TE, csMinusTrials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    RE.thirdOdor = extractReversalsFromTE(TE, Odor3Trials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    RE.csPlusReward = extractReversalsFromTE(TE, csPlusTrials & rewardTrials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    RE.allTrials = extractReversalsFromTE(TE, validTrials, dataToPull, 'maxReversals', 1000);%, 'maxReversals', 1);
    nReversals = size(RE.csPlus.phPeakPercentile_cs_ch1.after, 1);
    
if saveOn
    save(fullfile(savepath, ['RE_' subjectName '.mat']), 'RE');
    disp(['*** saving: ' fullfile(savepath, ['RE_' subjectName '.mat']) ' ***']);
end
    
%% reversal averages

smoothWindow = 5;

%     peakFieldCh1 = 'phPeakMean_cs_ch1';
%     peakFieldCh2 = 'phPeakMean_cs_ch2';    
    peakFieldCh1 = 'phPeakPercentile_cs_ch1';
    peakFieldCh2 = 'phPeakPercentile_cs_ch2';        
    saveName = [subjectName '_' strtok(strtok(peakFieldCh1, '_'), '_') '_phCue_revAvg'];  

    h=ensureFigure(saveName, 1);
    mcLandscapeFigSetup(h);
    subplot(2,2,1);

    % data fields from RE to plot
    dataFields = {
        'csLicks', 'k';...        
        peakFieldCh1, 'g';...
        };
    if ismember(2, channels)
        dataFields = [dataFields; {peakFieldCh2, 'r'}];
    end
    
    % CS MINUS -> CS PLUS
        % x data and indices for baseline (cs- baseline) and ceiling (cs+
        % baseline)
    xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
    bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);
    ceilIx(1) = nearest(RE.csPlus.trialsBefore, -30); ceilIx(2) = nearest(RE.csPlus.trialsBefore, 0);
    
    % scale data and plot
    nFields = size(dataFields, 1);
    for counter = 1:nFields
        dataField = dataFields{counter, 1};
        color = dataFields{counter, 2};
        revNorm = [RE.csMinus.(dataField).before RE.csPlus.(dataField).after];
        revNorm = smoothdata(revNorm, 2, 'movmean', smoothWindow, 'omitnan');
%         baseline = nanmean(RE.csMinus.(dataField).before(:,bl(1):bl(2)), 2); % subtract off Cs- prior to reversal
        ceiling = percentile(RE.csPlus.(dataField).before(:,ceilIx(1):ceilIx(2)), 0.8, 2); % divide by Cs+ prior to reversal
%         revNorm = revNorm - baseline; % RIGHT NOW JUST NORMALIZE, NO BASELINE SUBTRACTION    
        revNorm = revNorm ./ ceiling;
        subplot(2, 4, 1); hold on;
        boundedline(xData, nanmean(revNorm), nanSEM(revNorm), color);
        set(gca, 'XLim', [-50 50]);
        subplot(2,4,counter + 1);
        imagesc('XData', xData, 'CData', revNorm);
        set(gca, 'XLim', [-50 50]);
    end
%     set(gca, 'XLim', [-50 50]);
%     xlabel('Trials of new CS+ odor from reversal'); 
%     ylabel('Cue response normalized'); title([strtok(peakFieldCh1, '_') ', avg ' num2str(nReversals) ' reversals'], 'Interpreter', 'none');
    
    % CS PLUS -> CS MINUS
    subplot(2,1,2);
        % x data and indices for baseline (cs- baseline) and ceiling (cs+
        % baseline)        
    xData = [RE.csPlus.trialsBefore RE.csMinus.trialsAfter];
    bl(1) = nearest(RE.csMinus.trialsBefore, -30); bl(2) = nearest(RE.csMinus.trialsBefore, 0);
    ceilIx(1) = nearest(RE.csPlus.trialsBefore, -30); ceilIx(2) = nearest(RE.csPlus.trialsBefore, 0);
    
    % scale data and plot
    for counter = 1:size(dataFields, 1)
        dataField = dataFields{counter, 1};
        color = dataFields{counter, 2};
        revNorm = [RE.csPlus.(dataField).before RE.csMinus.(dataField).after];
        revNorm = smoothdata(revNorm, 2, 'movmean', smoothWindow, 'omitnan');
%         baseline = nanmean(RE.csMinus.(dataField).before(:,bl(1):bl(2)), 2); % subtract off Cs- prior to reversal
        ceiling = percentile(RE.csPlus.(dataField).before(:,ceilIx(1):ceilIx(2)), 0.8, 2); % divide by Cs+ prior to reversal
%         revNorm = revNorm - baseline; % RIGHT NOW JUST NORMALIZE, NO BASELINE SUBTRACTION         
        revNorm = revNorm ./ ceiling;
        subplot(2,4,5); hold on;
        boundedline(xData, nanmean(revNorm), nanSEM(revNorm), color);
        set(gca, 'XLim', [-50 50]);
        subplot(2,4,4 + counter + 1);
        imagesc('XData', xData, 'CData', revNorm);
        set(gca, 'XLim', [-50 50]);
    end
%     set(gca, 'XLim', [-50 50]);
%     xlabel('Trials of new CS- odor from reversal'); 
%     ylabel('Cue response normalized'); title([strtok(peakFieldCh1, '_') ', avg ' num2str(nReversals) ' reversals'], 'Interpreter', 'none');
    

    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
    
    