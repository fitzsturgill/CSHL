% 4/10/17  Analysis script for pavlovian reversals using LickNoLick_Odor_V2
% protocol



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
    BL{end + 1} = [0 4];
end
tau(1) = 1.5;
if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'expFit';
    BL{end + 1} = [0 4];    
end
tau(2) = 1;

%% baseline by trial
 TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL,...
    'tau', tau);   
%% baseline expfit
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', {'expFit', 'expFit'}, 'zeroField', 'Cue', 'channels', channels, 'baseline', BL,...
    'tau', tau);
% TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'expFitBegin', 0.1,...
%     'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL, 'downsample', 305);
%%
% if you are reloading TE do this:
channels = [1 2];

%% extract peak trial dFF responses to cues and reinforcement and lick counts
% csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);
csWindow = zeros(nTrials, 2);
csWindow(:,2) = cellfun(@(x,y,z) max(x(end), y(end)) - z(1), TE.AnswerLick, TE.AnswerNoLick, TE.Cue); 
% max 1) to select AnswerNoLick time stampfor no lick trials (unused state contains NaN)
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
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', 11, 'Fs', 20, 'startField', 'Start');

%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = uigetdir;
% basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\';
% basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_20and17_combined\';
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);


%%
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'trialType', 1));
%%
TE = addPupilometryToTE(TE, 'duration', 11, 'zeroField', 'Cue',  'frameRate', 60, 'frameRateNew', 20);
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
        figname2 = ['corrected_allTrials_ch1']; 
        ensureFigure(figname2, 1);
        [ha, hl] = phPlotAverageFromTE(TE, 1:length(TE.filename), 1,...
    'FluorDataField', 'ZS', 'window', [0.1, max(TE.Photometry.xData) - min(TE.Photometry.xData)], 'zeroTimes', TE.Photometry.startTime); %high value, reward
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
    Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
    Odor2Trials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);    
    
    rewardTrials = filterTE(TE, 'trialType', [1 5], 'reject', 0);
    hitTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    missTrials = filterTE(TE, 'trialOutcome', -1, 'reject', 0);
    punishTrials = filterTE(TE, 'trialType', [3 6], 'reject', 0);    
    neutralTrials = filterTE(TE, 'trialType', [2 4], 'reject', 0);
    block2Trials = filterTE(TE, 'BlockNumber', 2, 'reject', 0);
    block3Trials = filterTE(TE, 'BlockNumber', 3, 'reject', 0);
    csPlusTrials = filterTE(TE, 'trialType', [1 2], 'reject', 0);
    csMinusTrials = filterTE(TE, 'trialType', [3 4], 'reject', 0);
    uncuedReward = filterTE(TE, 'trialType', 5, 'reject', 0);
    uncuedPunish = filterTE(TE, 'trialType', 6, 'reject', 0);
    trialTypes = 1:6;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end

trialCount = [1:length(TE.filename)]';
%%
saveName = [subjectName '_longitudinalCueResponses_CS+'];
h=ensureFigure(saveName, 1); 
mcLandscapeFigSetup(h);
subplot(5,1,1); scatter(trialCount(csPlusTrials), TE.csLicks.rate(csPlusTrials),'o'); % both blLicks and anticipatoryLicks2 are from 2s periods
% subplot(4,1,1); scatter(1:length(allTE.epoch), cuedReward.anticipatoryLicks2);
maxLR = max(TE.csLicks.rate(csPlusTrials));
% subplot(4,1,1); plot(smooth(cuedReward.anticipatoryLicks1 / mean(cuedReward.blLicks)));
% maxLR = smooth(max(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxLR, 'g', 'Marker', 'none');
hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxLR, 'r', 'Marker', 'none');
ylabel('antic. licks'); title('CS+ and/or reward trials, dF/F is Zscored');
set(gca, 'XLim', [1 length(trialCount)]);

if ismember(1, channels)
    subplot(5,1,2); scatter(trialCount(csPlusTrials), TE.phPeakPercentile_cs(1).data(csPlusTrials), '.');
    maxP = max(TE.phPeakPercentile_cs(1).data(csPlusTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
    ylabel('BF: CS DF/F (50%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(2, channels)
    subplot(5,1,3); scatter(trialCount(csPlusTrials), TE.phPeakPercentile_cs(2).data(csPlusTrials), '.');
    maxP = max(TE.phPeakPercentile_cs(2).data(csPlusTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
ylabel('VTA:CS DF/F (50%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(1, channels)
    subplot(5,1,4); scatter(trialCount(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), '.');
    maxP = max(TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
    ylabel('BF: US DF/F (50%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(2, channels)
    subplot(5,1,5); scatter(trialCount(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials), '.');
    maxP = max(TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
ylabel('VTA: US DF/F (50%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
     xlabel('Trial Count');
end

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end


%% CS Minus
saveName = [subjectName '_longitudinalCueResponses_CSminus'];
h=ensureFigure(saveName, 1); 
mcLandscapeFigSetup(h);
subplot(5,1,1); scatter(trialCount(csMinusTrials), TE.csLicks.rate(csMinusTrials),'o'); % both blLicks and anticipatoryLicks2 are from 2s periods
% subplot(4,1,1); scatter(1:length(allTE.epoch), cuedReward.anticipatoryLicks2);
maxLR = max(TE.csLicks.rate(csMinusTrials));
% subplot(4,1,1); plot(smooth(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
% maxLR = smooth(max(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxLR, 'g', 'Marker', 'none');
hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxLR, 'r', 'Marker', 'none');
ylabel('antic. licks'); title('CS- and/or punish trials, dF/F is Zscored');
set(gca, 'XLim', [1 length(trialCount)]);

if ismember(1, channels)
    subplot(5,1,2); scatter(trialCount(csMinusTrials), TE.phPeakPercentile_cs(1).data(csMinusTrials), '.');
    maxP = max(TE.phPeakPercentile_cs(1).data(csMinusTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
    ylabel('BF: CS DF/F (50%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(2, channels)
    subplot(5,1,3); scatter(trialCount(csMinusTrials), TE.phPeakPercentile_cs(2).data(csMinusTrials), '.');
    maxP = max(TE.phPeakPercentile_cs(2).data(csMinusTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
ylabel('VTA:CS DF/F (50%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(1, channels)
    subplot(5,1,4); scatter(trialCount(csMinusTrials & punishTrials), TE.phPeakPercentile_us(1).data(csMinusTrials & punishTrials), '.');
    maxP = max(TE.phPeakPercentile_us(1).data(csMinusTrials & punishTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
    ylabel('BF: US DF/F (50%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(2, channels)
    subplot(5,1,5); scatter(trialCount(csMinusTrials & punishTrials), TE.phPeakPercentile_us(2).data(csMinusTrials & punishTrials), '.');
    maxP = max(TE.phPeakPercentile_us(2).data(csMinusTrials & punishTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
ylabel('VTA: US DF/F (50%)'); xlabel('Trial Count');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end
%% PH Rasters, CS+, CS-
CLimFactor = 2;
CLim = {[-0.005 0.005], [-0.1 0.1]};
CLim = {[-0.005 0.005], [-0.1 0.1]};
trialStart = TE.Photometry.xData(1);
reversals = find(TE.BlockChange);
for channel = channels

    saveName = [subjectName '_comboRasters_ch' num2str(channel)];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);
    

%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,4,1); 
    eventRasterFromTE(TE, csPlusTrials, 'Port1In', 'trialNumbering', 'global',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('CS+'); ylabel('trial number');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    

    
%     subplot(1,4,2); phRasterFromTE(TE, csPlusTrials, channel, 'CLimFactor', CLimFactor, 'CLim', CLim{channel}, 'trialNumbering', 'global');
%         set(gca, 'FontSize', 14)
    subplot(1,4,2); phRasterFromTE(TE, csPlusTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
        set(gca, 'FontSize', 14)
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines        
        
    subplot(1,4,3); 
    eventRasterFromTE(TE, csMinusTrials, 'Port1In', 'trialNumbering', 'global',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    title('CS-');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    
    subplot(1,4,4); phRasterFromTE(TE, csMinusTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    set(gca, 'FontSize', 14)
    xlabel('time from cue (s)');         
    colormap('parula');  

    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
end


%% photometry averages, zscored
%     ylim = [-2 8];
    saveName = [subjectName '_phAvgs'];  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);

    pm = [2 2];
    
    % - 6 0 4
    if ismember(1, channels)
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, neutralTrials, uncuedReward, uncuedPunish}, 1,...
            'FluorDataField', 'ZS', 'window', [3, 7], 'linespec', {'b', 'r', 'k', 'c', 'm'}); %high value, reward
        legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);
    end
    
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, neutralTrials, uncuedReward, uncuedPunish}, 2,...
            'FluorDataField', 'ZS', 'window', [3, 7], 'linespec', {'b', 'r', 'k', 'c', 'm'}); %high value, reward
        legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        ylabel('VTA dF/F Zscored'); xlabel('time from cue (s)'); %set(gca, 'YLim', ylim);
    end
    
    % - 6 0 4
    if ismember(1, channels)    
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 1,...
        'FluorDataField', 'ZS', 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('CS+, outcomes'); set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);
    end
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 2,...
            'FluorDataField', 'ZS', 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        xlabel('time from cue (s)');     set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);
    end
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end    
    
    %%
    saveName = [subjectName '_scatter'];
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);
    pm = [2 2];    
    if all(ismember([1 2], channels))    
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1);
        scatter(TE.phPeakPercentile_cs(2).data(csPlusTrials), TE.phPeakPercentile_cs(1).data(csPlusTrials), '.'); hold on;
        fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
        fob = fit(TE.phPeakPercentile_cs(2).data(csPlusTrials), TE.phPeakPercentile_cs(1).data(csPlusTrials), 'poly1', fo); 
        plot(fob,'predfunc'); legend off;
        title('CS+, Cue'); ylabel('BF dF/F Zscored'); xlabel('VTA dF/F'); textBox(subjectName);
    end

    if all(ismember([1 2], channels))    
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1);
        scatter(TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), '.'); hold on;
        fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
        fob = fit(TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), 'poly1', fo); 
        plot(fob,'predfunc'); legend off;    
        title('CS+, Reward'); ylabel('BF dF/F Zscored'); xlabel('VTA dF/F');  
    end
    
    if ismember(1, channels)    
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1);
        scatter(TE.csLicks.rate(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), 'o'); hold on;
        fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
        try
            fob = fit(TE.csLicks.rate(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), 'poly1', fo); 
            plot(fob,'predfunc'); legend off;      
        end
        title('CS+, reward vs. cue licks'); ylabel('BF dF/F Zscored'); xlabel('Antic. Lick Rate');    
    end
    
    if ismember(2, channels)        
        subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1);
        scatter(TE.csLicks.rate(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials), 'o'); hold on;
        fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
        try
            fob = fit(TE.csLicks.rate(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials), 'poly1', fo); 
            plot(fob,'predfunc'); legend off;         
        end
        ylabel('VTA dF/F Zscored'); xlabel('Antic. Lick Rate');
    end
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end    

 %% single trial traces
    nTraces = 20;
    saveName = [subjectName '_singleTrial3'];  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);
    csPlusRewardTrials = find(csPlusTrials & rewardTrials & hitTrials);
    rn = rand(length(csPlusRewardTrials), 1);
    [~, I] = sort(rn);
    %
    subplot(2,1,1, 'FontSize', 12, 'LineWidth', 1); plot(TE.Photometry.xData, TE.Photometry.data(1).dFF(I(1:nTraces), :)', 'k'); hold on;
    phPlotAverageFromTE(TE, csPlusTrials & rewardTrials & hitTrials, 1, 'window', [-4, 7], 'linespec', {'r'});
    set(gca, 'XLim', [-4, 7]); ylabel('BF dF/F'); xlabel('time from cue (s)');
    
    
    subplot(2,1,2, 'FontSize', 12, 'LineWidth', 1); plot(TE.Photometry.xData, TE.Photometry.data(1).dFF(I(1:5), :)', 'k');
    set(gca, 'XLim', [-4, 7]); ylabel('BF dF/F'); xlabel('time from cue (s)');
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end    
    
    %% post-reversal analysis, assuming you always have photometry in channel 1, but not necessarily in channel 2
    
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
        'trialNumber', TE.trialNumber...
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
    RE.csPlus = extractReversalsFromTE(TE, csPlusTrials, dataToPull, 'maxReversals', 1);
    RE.csMinus = extractReversalsFromTE(TE, csMinusTrials, dataToPull, 'maxReversals', 1);
    RE.csPlusReward = extractReversalsFromTE(TE, csPlusTrials & rewardTrials, dataToPull, 'maxReversals', 1);    
    nReversals = size(RE.csPlus.phPeakPercentile_cs_ch1.after, 1);
    
if saveOn
    save(fullfile(savepath, ['RE_' subjectName '.mat']), 'RE');
    disp(['*** saving: ' fullfile(savepath, ['RE_' subjectName '.mat']) ' ***']);
end
    
    %% reversal averages

    peakFieldCh1 = 'phPeakMean_cs_ch1';
    peakFieldCh2 = 'phPeakMean_cs_ch2';    
%     peakFieldCh1 = 'phPeakPercentile_cs_ch1';
%     peakFieldCh2 = 'phPeakPercentile_cs_ch2';        
    saveName = [subjectName '_' strtok(strtok(peakFieldCh1, '_'), '_') '_phCue_revAvg'];  
    h=ensureFigure(saveName, 1);
    mcLandscapeFigSetup(h);
    subplot(2,1,1);
    xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
    bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);bl(3) = nearest(xData, 20);
    revNormAvg_ch1 = [RE.csMinus.(peakFieldCh1).before RE.csPlus.(peakFieldCh1).after];
%     revNormAvg_ch1 = nanfastsmooth(nanmean(revNormAvg_ch1), 5);
    revNormAvg_ch1 = nanmean(revNormAvg_ch1);
%     revNormAvg_ch1 = nanmean(revNormAvg_ch1);
%     revNormAvg_ch1 = revNormAvg_ch1 - nanmean(revNormAvg_ch1(bl(1):bl(2)));
%     revNormAvg_ch1 = revNormAvg_ch1 / percentile(revNormAvg_ch1(bl(2):bl(3)), 0.90);
    
    revNormAvg_ch2 = [RE.csMinus.(peakFieldCh2).before RE.csPlus.(peakFieldCh2).after];
%     revNormAvg_ch2 = nanfastsmooth(nanmean(revNormAvg_ch2), 5);
    revNormAvg_ch2 = nanmean(revNormAvg_ch2);    
%     revNormAvg_ch2 = nanmean(revNormAvg_ch2);
%     revNormAvg_ch2 = revNormAvg_ch2 - nanmean(revNormAvg_ch2(bl(1):bl(2)));
%     revNormAvg_ch2 = revNormAvg_ch2 / percentile(revNormAvg_ch2(bl(2):bl(3)), 0.90);

    plot(xData, smooth(revNormAvg_ch1), 'g'); hold on;
    plot(xData, smooth(revNormAvg_ch2), 'r');
    set(gca, 'XLim', [-10 20], 'YLim', [-2 2]);xlabel('Trials of new CS+ odor from reversal'); 
    ylabel('Cue dFF ZScored'); title([strtok(peakFieldCh1, '_') ', avg ' num2str(nReversals) ' reversals'], 'Interpreter', 'none');

    % reinforcment response
    subplot(2,1,2);
    xData = [RE.csPlusReward.trialsBefore RE.csPlusReward.trialsAfter];
    bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);bl(3) = nearest(xData, 80);
    revRewNormAvg_ch1 = [RE.csPlusReward.phPeakPercentile_us_ch1.before RE.csPlusReward.phPeakPercentile_us_ch1.after];
%     revRewNormAvg_ch1 = nanfastsmooth(nanmean(revRewNormAvg_ch1), 5);
    revRewNormAvg_ch1 = nanmean(revRewNormAvg_ch1);
%     revRewNormAvg_ch1 = revRewNormAvg_ch1 - nanmean(revRewNormAvg_ch1(bl(1):bl(2)));
%     revRewNormAvg_ch1 = revRewNormAvg_ch1 / percentile(revRewNormAvg_ch1(bl(2):bl(3)), 0.90);
    
    revRewNormAvg_ch2 = [RE.csPlusReward.(peakFieldCh2).before RE.csPlusReward.(peakFieldCh2).after];
%     revRewNormAvg_ch2 = nanfastsmooth(nanmean(revRewNormAvg_ch2), 5);
    revRewNormAvg_ch2 = nanmean(revRewNormAvg_ch2);
%     revRewNormAvg_ch2 = revRewNormAvg_ch2 - nanmean(revRewNormAvg_ch2(bl(1):bl(2)));
%     revRewNormAvg_ch2 = revRewNormAvg_ch2 / max(revRewNormAvg_ch2); %percentile(revRewNormAvg_ch2(bl(2):bl(3)), 0.90);

%     revRewNormAvg_ch2 = revRewNormAvg_ch2 / percentile(revRewNormAvg_ch2(bl(2):bl(3)), 0.90);

    plot(xData, revRewNormAvg_ch1, 'g'); hold on;
    plot(xData, revRewNormAvg_ch2, 'r');
    set(gca, 'XLim', [-40 80]);xlabel('Trials of CS+ odor from reversal');
    ylabel('Reward dFF ZScored'); title([strtok(peakFieldCh1, '_') ', avg ' num2str(nReversals) ' reversals'], 'Interpreter', 'none');
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
    
    
    %% phCue Reversal array
    saveName = [subjectName '_' strtok(strtok(peakFieldCh1, '_'), '_') '_phCue_revArray'];  
    h = ensureFigure(saveName, 1);
    mcLandscapeFigSetup(h);
    nReversals = size(RE.csPlus.(peakFieldCh1).after, 1);
    ass = ceil(sqrt(nReversals));
    xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
    bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);bl(3) = nearest(xData, 0);
    for counter = 1:nReversals
        subplot(ass,ass,counter);
        revNorm_ch1 = [RE.csMinus.(peakFieldCh1).before(counter,:) RE.csPlus.(peakFieldCh1).after(counter,:)];
        revNorm_ch1 = nanfastsmooth(revNorm_ch1, 5,1);
%         revNorm_ch1 = revNorm_ch1 - nanmean(revNorm_ch1(bl(1):bl(2)));
%         revNorm_ch1 = revNorm_ch1 / percentile(revNorm_ch1(bl(2):end), 0.90);
        rev_phBaseline_ch1 = [RE.csMinus.phBaseline_ch1.before(counter,:) RE.csPlus.phBaseline_ch1.after(counter,:)];   
        rev_phBaseline_ch1 = nanfastsmooth(rev_phBaseline_ch1, 5, 1);

    
        revNorm_ch2 = [RE.csMinus.(peakFieldCh2).before(counter,:) RE.csPlus.(peakFieldCh2).after(counter,:)];
        revNorm_ch2 = nanfastsmooth(revNorm_ch2, 5,1);
%         revNorm_ch2 = revNorm_ch2 - nanmean(revNorm_ch2(bl(1):bl(2)));
%         revNorm_ch2 = revNorm_ch2 / percentile(revNorm_ch2(bl(2):end), 0.90);
        rev_phBaseline_ch2 = [RE.csMinus.phBaseline_ch2.before(counter,:) RE.csPlus.phBaseline_ch2.after(counter,:)];   
        rev_phBaseline_ch2 = nanfastsmooth(rev_phBaseline_ch2, 5, 1);

        revCsLicks = [RE.csMinus.csLicks.before(counter,:) RE.csPlus.csLicks.after(counter,:)];
%         revCsLicks = nanfastsmooth(revCsLicks, 5,1);

%         plot(xData, revNorm_ch1, 'b'); hold on;
%         plot(xData, revNorm_ch2, 'r'); hold on;
        plot(xData, rev_phBaseline_ch1, 'c'); hold on;
        plot(xData, rev_phBaseline_ch2, 'm');        
%         plot(xData, revCsLicks, 'b');
        set(gca, 'XLim', [-40 80]); 
        if counter == 1
            title(['Peak method: ' peakFieldCh1], 'Interpreter', 'none');
            xlabel('Trials of new CS+ odor from reversal');
            ylabel('Cue dF Zscored'); 
        end
    end
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
    
      %% phOutcome Reversal array
    saveName = [subjectName '_' strtok(strtok(peakFieldCh1, '_'), '_') '_phOutcome_revArray'];  
    h = ensureFigure(saveName, 1);
    mcLandscapeFigSetup(h);
    nReversals = size(RE.csPlus.(peakFieldCh1).after, 1);
    ass = ceil(sqrt(nReversals));
    xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
    bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);bl(3) = nearest(xData, 0);
    for counter = 1:nReversals
        subplot(ass,ass,counter);
        revRewNormAvg_ch1 = [RE.csPlusReward.phPeakPercentile_us_ch1.before RE.csPlusReward.phPeakPercentile_us_ch1.after];
        
        revNorm_ch1 = [RE.csMinus.(peakFieldCh1).before(counter,:) RE.csPlus.(peakFieldCh1).after(counter,:)];
        revNorm_ch1 = nanfastsmooth(revNorm_ch1, 5,1);
%         revNorm_ch1 = revNorm_ch1 - nanmean(revNorm_ch1(bl(1):bl(2)));
%         revNorm_ch1 = revNorm_ch1 / percentile(revNorm_ch1(bl(2):end), 0.90);
        rev_phBaseline_ch1 = [RE.csMinus.phBaseline_ch1.before(counter,:) RE.csPlus.phBaseline_ch1.after(counter,:)];   
        rev_phBaseline_ch1 = nanfastsmooth(rev_phBaseline_ch1, 5, 1);

    
        revNorm_ch2 = [RE.csMinus.(peakFieldCh2).before(counter,:) RE.csPlus.(peakFieldCh2).after(counter,:)];
        revNorm_ch2 = nanfastsmooth(revNorm_ch2, 5,1);
%         revNorm_ch2 = revNorm_ch2 - nanmean(revNorm_ch2(bl(1):bl(2)));
%         revNorm_ch2 = revNorm_ch2 / percentile(revNorm_ch2(bl(2):end), 0.90);
        rev_phBaseline_ch2 = [RE.csMinus.phBaseline_ch2.before(counter,:) RE.csPlus.phBaseline_ch2.after(counter,:)];   
        rev_phBaseline_ch2 = nanfastsmooth(rev_phBaseline_ch2, 5, 1);

        revCsLicks = [RE.csMinus.csLicks.before(counter,:) RE.csPlus.csLicks.after(counter,:)];
%         revCsLicks = nanfastsmooth(revCsLicks, 5,1);

%         plot(xData, revNorm_ch1, 'b'); hold on;
%         plot(xData, revNorm_ch2, 'r'); hold on;
        plot(xData, rev_phBaseline_ch1, 'c'); hold on;
        plot(xData, rev_phBaseline_ch2, 'm');        
%         plot(xData, revCsLicks, 'b');
        set(gca, 'XLim', [-40 80]); 
        if counter == 1
            title(['Peak method: ' peakFieldCh1], 'Interpreter', 'none');
            xlabel('Trials of new CS+ odor from reversal');
            ylabel('Cue dF Zscored'); 
        end
    end
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
    
    %% Cue Lick Reversal Array
    saveName = [subjectName '_cueLicks_revArray'];  
    h = ensureFigure(saveName, 1);
    mcLandscapeFigSetup(h);
    nReversals = size(RE.csPlus.csLicks.after, 1);
    ass = ceil(sqrt(nReversals));
    xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
    bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);bl(3) = nearest(xData, 0);
    for counter = 1:nReversals
        subplot(ass,ass,counter);
        revCsLicks = [RE.csMinus.csLicks.before(counter,:) RE.csPlus.csLicks.after(counter,:)];
        revCsLicks = nanfastsmooth(revCsLicks, 5,1);
        
%         revNorm_ch2 = revNorm_ch2 - nanmean(revNorm_ch2(bl(1):bl(2)));
%         revNorm_ch2 = revNorm_ch2 / percentile(revNorm_ch2(bl(2):end), 0.90);

        plot(xData, revCsLicks, 'g'); hold on;

        set(gca, 'XLim', [-40 80]); 
        if counter == 1
            title('CS Licks', 'Interpreter', 'none');
            xlabel('Trials of new CS+ odor from reversal');
            ylabel('Licks/s'); 
        end
    end
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
    
    %% trial - 1 vs trial correlation (diff)
    ensureFigure('trial_diff', 1);
    ch1Data = diff([RE.csMinus.(peakFieldCh1).before RE.csPlus.(peakFieldCh1).after], 2);
    ch2Data = diff([RE.csMinus.(peakFieldCh2).before RE.csPlus.(peakFieldCh2).after], 2);
    ch1Data = reshape(ch1Data, numel(ch1Data), 1);
    ch2Data = reshape(ch2Data, numel(ch2Data), 1);
    ch1Data = ch1Data(~isnan(ch1Data));
    ch2Data = ch2Data(~isnan(ch2Data));
    
    scatter(ch2Data, ch1Data, '.'); hold on;
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(ch2Data, ch1Data, 'poly1', fo); 
    plot(fob,'predfunc'); legend off;
    addUnityLine;
    title('CS+, Cue, trial = n vs n + 1 difference'); ylabel('BF delta zscore'); xlabel('VTA delta zscore'); textBox(subjectName);
    
    
   %% accuracy split
    saveName = [subjectName '_' strtok(strtok(peakFieldCh1, '_'), '_') '_phCue_revAvg_accuracySplit'];  
    h=ensureFigure(saveName, 1);
    mcLandscapeFigSetup(h);
    xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
    accuracy_before = sum(ismember(RE.csMinus.trialOutcome.before, [1 2]), 2) ./ sum(~isnan(RE.csMinus.trialOutcome.before), 2);
    [~, I] = sort(accuracy_before);   
    bottom = I(1:floor(numel(I) / 2));
    top = I(floor(numel(I) / 2) + 1:end);

    bottom_revNormAvg_ch1= nanfastsmooth(nanmean([RE.csMinus.(peakFieldCh1).before(bottom,:) RE.csPlus.(peakFieldCh1).after(bottom,:)]), 3);
    top_revNormAvg_ch1= nanfastsmooth(nanmean([RE.csMinus.(peakFieldCh1).before(top,:) RE.csPlus.(peakFieldCh1).after(top,:)]), 3);

    bottom_revNormAvg_ch2= nanfastsmooth(nanmean([RE.csMinus.(peakFieldCh2).before(bottom,:) RE.csPlus.(peakFieldCh2).after(bottom,:)]), 3);
    top_revNormAvg_ch2= nanfastsmooth(nanmean([RE.csMinus.(peakFieldCh2).before(top,:) RE.csPlus.(peakFieldCh2).after(top,:)]), 3);
    
    subplot(2,2,1);
    plot(xData, bottom_revNormAvg_ch1, 'g--'); hold on;
    plot(xData, top_revNormAvg_ch1, 'g');    
    plot(xData, bottom_revNormAvg_ch2, 'r--'); hold on;
    plot(xData, top_revNormAvg_ch2, 'r');    
    legend({'bottom', 'top', 'bottom', 'top'}, 'Location', 'southeast', 'Box', 'off');
    set(gca, 'XLim', [-40 80]);xlabel('Trials of CS+ odor from reversal');
    ylabel('Cue dFF ZScored'); title([strtok(peakFieldCh1, '_') ', avg ' num2str(nReversals) ' reversals'], 'Interpreter', 'none');
    subplot(2,2,2);
    plot(xData, bottom_revNormAvg_ch1, 'g--'); hold on;
    plot(xData, top_revNormAvg_ch1, 'g'); 
    legend({'bottom', 'top'}, 'Location', 'southeast', 'Box', 'off');
    set(gca, 'XLim', [-40 80]);xlabel('Trials of CS+ odor from reversal');
    ylabel('Cue dFF ZScored, BF')
    title('Bottom vs top 50% of pre-reversal acurracy');
    subplot(2,2,3);
    plot(xData, bottom_revNormAvg_ch2, 'r--'); hold on;
    plot(xData, top_revNormAvg_ch2, 'r'); 
    legend({'bottom', 'top'}, 'Location', 'southeast', 'Box', 'off');
    set(gca, 'XLim', [-40 80]);xlabel('Trials of CS+ odor from reversal');
    ylabel('Cue dFF ZScored, VTA')

    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
    
    %%
    ensureFigure('test', 1);
    subplot(1,3,1);
    title('Velocity');
    image(TE.Wheel.data.V(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
%     subplot(3,2,1); image(rewards_chat, 'XData', window, 'CDataMapping', 'Scaled'); 
%     set(gca, 'CLim', [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); 
    colormap('jet');  
    title('Velocity');    
    subplot(1,3,2); phRasterFromTE(TE, csPlusTrials, 1, 'trialNumbering', 'consecutive'); % 'CLimFactor', CLimFactor,
    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,3,3); phRasterFromTE(TE, csPlusTrials, 2, 'trialNumbering', 'consecutive'); % 'CLimFactor', CLimFactor,
    title('DAT');    
    
    
    
    
    