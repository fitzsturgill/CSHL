% uncertainty blocks analysis script suing LickNoLick_Odor_V2 protocol

saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_LNL_odor_V2(sessions);
%%
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRGECO1a (ch2)
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
    BL{end + 1} = [0 4];    
end


    

% TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'expFit', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL);
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'expFitBegin', 0.1,...
    'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL, 'downsample', 10);

%%
% if you are reloading TE do this:
channels = [1 2];

%% extract peak trial dFF responses to cues and reinforcement and lick counts
% csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);
csWindow = zeros(nTrials, 2);
csWindow(:,2) = cellfun(@(x,y,z) max(x(end), y(end)) - z(1), TE.AnswerLick, TE.AnswerNoLick, TE.Cue); 
% max 1) to select AnswerNoLick time stamp for no lick trials (unused state contains NaN)
% 2) to select AnswerLick time stamp for lick trials (AnswerLick follows
% AnswerNoLick state)

ch1CsWindow = [0.25 1];
ch2CsWindow = [0.25 1];

TE.csLicks = countEventFromTE(TE, 'Port1In', csWindow, TE.Cue);

usWindow = [0 0.5];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'


for channel = channels
    TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'mean', 'phField', 'ZS');
    if channel == 1
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');
    elseif channel == 2
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');        
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');  
    end
end



TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);
%%
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'trialType', 1));

%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = uigetdir;
% basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\';
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

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
    if channel == 1
        figname = ['trialBleach_Correction_ch' num2str(channel)];
        ensureFigure(figname, 1);
        nSessions = length(TE.Photometry.bleachFit);
        subA = ceil(sqrt(nSessions));
        for counter = 1:nSessions
            subplot(subA, subA, counter);
            plot(TE.Photometry.bleachFit(counter, channel).trialTemplateFullX, TE.Photometry.bleachFit(counter, channel).trialTemplateFull, 'g'); hold on;            
            plot(TE.Photometry.bleachFit(counter, channel).trialTemplateX, TE.Photometry.bleachFit(counter, channel).trialTemplate, 'b'); 
            plot(TE.Photometry.bleachFit(counter, channel).fitX, TE.Photometry.bleachFit(counter, channel).trialFit, 'r');
        %     title(num2str(counter));    
        end
        % average of all trials for this channel to eyeball correction
        figname2 = ['corrected_allTrials_ch1']; 
        ensureFigure(figname2, 1);
        subplot(1,2,1);
        [ha, hl] = phPlotAverageFromTE(TE, 1:length(TE.filename), 1,...
    'FluorDataField', 'ZS', 'window', [0.1, max(TE.Photometry.xData) - min(TE.Photometry.xData)], 'zeroTimes', TE.Photometry.startTime); 
        title('corrected'); 
        set(gca, 'XLim', [0 9]); % kludge
        subplot(1,2,2);
        [ha, hl] = phPlotAverageFromTE(TE, 1:length(TE.filename), 1,...
    'FluorDataField', 'raw', 'window', [0.1, max(TE.Photometry.xData) - min(TE.Photometry.xData)], 'zeroTimes', TE.Photometry.startTime); 
        title('raw'); 
        set(gca, 'XLim', [0 9]); % kludge
        if saveOn
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.fig'));
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.jpg'));
            saveas(gcf, fullfile(savepath, [figname2 '.fig']));
            saveas(gcf, fullfile(savepath, [figname2 '.jpg']));            
        end
    end
end

%% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks    blocks 2 and 3
    validTrials = filterTE(TE, 'reject', 0);
    highTrials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
    mediumTrials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);    
    lowTrials = filterTE(TE, 'OdorValveIndex', 3, 'reject', 0);    
    rewardTrials = filterTE(TE, 'trialType', [1 3 5 7], 'reject', 0);
    neutralTrials = filterTE(TE, 'trialType', [2 4 6], 'reject', 0);
    trialTypes = 1:7;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end
    trialCount = (1:length(TE.filename))';
    
%% Rasters by cue condition
CLimFactor = 3;
CLim = {[-0.005 0.005], [-0.1 0.1]};
CLim = {[-0.005 0.005], [-0.1 0.1]};
trialStart = TE.Photometry.xData(1);
for channel = channels

    saveName = [subjectName '_cueRasters_ch' num2str(channel)];
    h=ensureFigure(saveName, 1);
    mcLandscapeFigSetup(h);
    
    % pReward = high
    % licks
    subplot(1,6,1); 
    eventRasterFromTE(TE, highTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('pRew = 0.75'); ylabel('trial number');
    set(gca, 'XLim', [-4 3]); 
%     set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    % photometry
    subplot(1,6,2); phRasterFromTE(TE, highTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-4 3]);
        set(gca, 'FontSize', 14);     xlabel('time from cue (s)');        

    % pReward = medium
    % licks
    subplot(1,6,3); 
    eventRasterFromTE(TE, mediumTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('pRew = 0.5');
    set(gca, 'XLim', [-4 3]); 
%     set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    % photometry
    subplot(1,6,4); phRasterFromTE(TE, mediumTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-4 3]);
        set(gca, 'FontSize', 14)
        
    % pReward = low        
    % licks
    subplot(1,6,5); 
    eventRasterFromTE(TE, lowTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('pRew = 0.25'); 
    set(gca, 'XLim', [-4 3]); 
%     set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    % photometry
    subplot(1,6,6); phRasterFromTE(TE, lowTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'window', [-4 3]);
    set(gca, 'FontSize', 14)        
 
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
end

%% Photometry averages- just do for channel 1 for now...
    h=ensureFigure('Photometry_Averages', 1); 
    mcLandscapeFigSetup(h);

    pm = [3 2];
    
    % First row: by cue condition
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-3 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
    axh = [];
    subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, {highTrials, mediumTrials, lowTrials}, 'Port1In', varargin{:});
    legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    title('Licking'); ylabel('Cue only'); textBox(TE.filename{1}(1:6)); set(gca, 'XLim', [-3 3]);
    
    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, {highTrials, mediumTrials, lowTrials}, 1, 'FluorDataField', 'ZS', 'window', [-3 3]); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    title('Photometry'); set(gca, 'XLim', [-3 3]);
    
    % Second row: Reward trials
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-3 7], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
    axh = [];
    subplot(pm(1), pm(2), 3); [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 3 5 7]), 'Port1In', varargin{:});
    legend(hl, {'pRhigh', 'pRmedium', 'pRlow', 'uncued'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    ylabel('Reward trials'); set(gca, 'XLim', [-3 7]);
    
    subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 3 5 7]), 1, 'FluorDataField', 'ZS', 'window', [-3 7]); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow', 'uncued'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    set(gca, 'XLim', [-3 7]);
    % Third row: Neutral trials
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-3 7], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
    axh = [];
    subplot(pm(1), pm(2), 5); [ha, hl] = plotEventAverageFromTE(TE, trialsByType([2 4 6]), 'Port1In', varargin{:});
    legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    ylabel('Neutral trials'); xlabel('time from cue (s)');set(gca, 'XLim', [-3 7]);
    
    subplot(pm(1), pm(2), 6, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([2 4 6]), 1, 'FluorDataField', 'ZS', 'window', [-3 7]); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    xlabel('time from cue (s)');set(gca, 'XLim', [-3 7]);
    
if saveOn    
    saveas(gcf, fullfile(savepath, 'phAverages.fig'));
    saveas(gcf, fullfile(savepath, 'phAverages.jpg'));
end
