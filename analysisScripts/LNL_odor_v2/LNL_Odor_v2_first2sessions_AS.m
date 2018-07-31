% LNL_Odor_v2_first2sessions_AS
% 4/10/17  Analysis script for pavlovian reversals using LickNoLick_Odor_V2
% protocol

saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_LNL_odor_V2(sessions);
%%
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRCaMP or jRGECO (ch2)
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
 TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL,...
    'tau', 2);   

% if you are reloading TE do this:
channels = [1 2];
%% extract cue and outcome responses, allowing for potentially increased trace period across 1st and 2nd sessions

% end of cue -> end of trace period!!!
csWindow = [cellfun(@(x) x(end), TE.Cue) cellfun(@(y,z) max(y(end), z(end)), TE.AnswerLick, TE.AnswerNoLick)]; % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);


usWindow = [0 0.75];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'

% percentile value for peak estimations
percentValue = 0.8;

% estimate respones for different events in each trial for photometry
% (bpCalcPeak_dFF) and for licking (countEventFromTE)
TE.csLicks = countEventFromTE(TE, 'Port1In', csWindow);
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);



for channel = channels
    TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'mean', 'phField', 'ZS');
    if channel == 1
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow - TE.Photometry.startTime, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow - TE.Photometry.startTime, TE.Photometry.startTime, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');
    elseif channel == 2
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow - TE.Photometry.startTime, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow - TE.Photometry.startTime, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');        
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'percentile', 'percentile', percentValue, 'phField', 'ZS');  
    end
end

%%
basepath = uigetdir; % prompts windows/mac osx to give you a location to save
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
% %     figname = ['sessionBleach_Correction_ch' num2str(channel)];
% %     ensureFigure(figname, 1);
% %     plot(TE.Photometry.data(channel).blF_raw, 'k'); hold on;
% %     plot(TE.Photometry.data(channel).blF, 'r');
% %     if saveOn
% %         saveas(gcf, fullfile(savepath, [figname '.fig']));
% %         saveas(gcf, fullfile(savepath, [figname '.jpg']));
% %     end
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
% find first reversal if it occurs on day 2
firstReversalBlock = TE.BlockNumber(find(TE.sessionIndex == 2 & TE.OdorValveIndex == 1 & TE.CSValence == -1, 1)); % block number where 1st odor has a negative valence.
firstReversalTrial = find(TE.BlockNumber == firstReversalBlock, 1);
day1Trials = TE.sessionIndex == 1 & validTrials;
day2Trials = TE.sessionIndex == 2 & validTrials;
day2Trials(firstReversalTrial:end) = false; % if firstReversalTrial is empty nothing happens...

%% photometry and lick averages- compare day 1 and day 2
%     ylim = [-2 8];
    saveName = [subjectName '_phAvgs'];  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);

    pm = [2 2];    
    window = [-4 3];
    % photometry averages
    if ismember(1, channels)
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & day1Trials & rewardTrials, csPlusTrials & day2Trials & rewardTrials}, 1,...
            'FluorDataField', 'ZS', 'window', window, 'linespec', {'c','b'}, 'zeroTimes', TE.Cue); %high value, reward
        legend(hl, {'cued reward, day 1', 'cued reward, day 2'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Acquisition'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);
    end
    
    if ismember(2, channels)
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & day1Trials & rewardTrials, csPlusTrials & day2Trials & rewardTrials}, 2,...
            'FluorDataField', 'ZS', 'window', window, 'linespec', {'c','b'}, 'zeroTimes', TE.Cue); %high value, reward
%         legend(hl, {'cuedReward, day 1', 'cuedReward, day 2'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        ylabel('VTA dF/F Zscored'); xlabel('time from cue (s)'); textBox(subjectName);%set(gca, 'YLim', ylim); 
    end
    
    % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', window, 'zeroTimes', TE.Cue,'linespec', {'c','b'}};
    axh = [];
    subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = plotEventAverageFromTE(TE, {csPlusTrials & day1Trials & rewardTrials, csPlusTrials & day2Trials & rewardTrials}, 'Port1In', varargin{:});
    ylabel('licks (s)'); xlabel('time from cue (s)'); 
    
    
%% Raster Plots for licking and both photometry channels
% CLimFactor = 2;
% 
% subplot(1,4,1); 
% eventRasterFromTE(TE, cuedReward, 'Port1In', 'trialNumbering', 'consecutive',...
%     'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
% title('cued'); ylabel('trial number');
% set(gca, 'YLim', [0 sum(cuedReward)]);
% set(gca, 'XLim', [-4 7]); 
% set(gca, 'FontSize', 14)
% 
% 
% saveName = 'Rasters_first2Sessions';
% h=ensureFigure(saveName, 1);
% mcPortraitFigSetup(h);
%     
% subplot(1,4,2);
% phRasterFromTE(TE, cuedReward, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
% set(gca, 'FontSize', 14);
%    
% 
%     
%     subplot(1,4,2); phRasterFromTE(TE, cuedReward, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
%         set(gca, 'FontSize', 14);
% 
%     subplot(1,4,3); phRasterFromTE(TE, omit, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
%         set(gca, 'FontSize', 14); title('omission');
% 
%     subplot(1,4,4); phRasterFromTE(TE, uncuedReward, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
%         set(gca, 'FontSize', 14); title('uncued');       
% 
%     if saveOn
%         saveas(gcf, fullfile(savepath, [saveName '.fig']));
%         saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
%     end
% end

