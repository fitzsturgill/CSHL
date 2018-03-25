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
    BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
    BL{end + 1} = [0 4];    
end

%% baseline by trial
 TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL,...
    'tau', 2);   
%% baseline expfit
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', {'expFit', 'expFit'}, 'zeroField', 'Cue', 'channels', channels, 'baseline', BL,...
    'tau', 2);
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

TE.csLicks = countEventFromTE(TE, 'Port1In', [0 3], TE.Cue);

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
pupLag = 0.3;
TE.wheelBaseline = mean(TE.Wheel.data.V(:,bpX2pnt(-3, 20, -4):bpX2pnt(0, 20, -4)), 2);
TE.pupilBaseline = mean(TE.pupil.pupDiameterNorm(:,bpX2pnt(-3, 20, -4):bpX2pnt(0, 20, -4)), 2);
% ch1CsWindow = [0.25 1];
TE.pupil_cs = mean(TE.pupil.pupDiameterNorm(:,bpX2pnt(ch1CsWindow(1) + pupLag, 20, -4):bpX2pnt(ch1CsWindow(2) + pupLag, 20, -4)), 2);
% TE.wheel_cs = mean(TE.Wheel.data.V(:,bpX2pnt(ch1CsWindow(1) + pupLag, 20, -4):bpX2pnt(ch1CsWindow(2) + pupLag, 20, -4)), 2);


TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);

%%
phasicWindow = [-2.7, -2.2];
sustainedWindow = [-1.5, 0];
TE.phasicAvg= bpCalcPeak_dFF(TE.Photometry, 1, phasicWindow, TE.Us, 'method', 'mean');
TE.sustainedAvg = bpCalcPeak_dFF(TE.Photometry, 1, sustainedWindow, TE.Us, 'method', 'mean');
TE.phasicLicks = countEventFromTE(TE, 'Port1In', phasicWindow, TE.Us);
TE.sustainedLicks = countEventFromTE(TE, 'Port1In', sustainedWindow, TE.Us);

pupLag = 0.3; % pupil dilation lags cholinergic signal (Nelson and Mooney)
TE.phasicPup = bpCalcPeak_Pupil(TE.pupil, phasicWindow + pupLag, TE.Us);
TE.sustainedPup = bpCalcPeak_Pupil(TE.pupil, sustainedWindow + pupLag, TE.Us);
TE.pupBaseline = bpCalcPeak_Pupil(TE.pupil, [1 4]);

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
    if channel == 1
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
        set(gca, 'XLim', [3 6.9]);
    end
    
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, neutralTrials, uncuedReward, uncuedPunish}, 2,...
            'FluorDataField', 'ZS', 'window', [3, 7], 'linespec', {'b', 'r', 'k', 'c', 'm'}); %high value, reward
        legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        ylabel('VTA dF/F Zscored'); xlabel('time from cue (s)'); %set(gca, 'YLim', ylim);
        set(gca, 'XLim', [3 6.9]);
    end
    
    lickThresh = 1;
    expectTrials = TE.csLicks.rate > lickThresh;
    surpriseTrials = TE.csLicks.rate <= lickThresh;
    % - 6 0 4
    if ismember(1, channels)    
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & expectTrials, csPlusTrials & rewardTrials & surpriseTrials}, 1,...
        'FluorDataField', 'ZS', 'window', [-3, 6.9], 'linespec', {'c', 'm'}); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('CS+, outcomes'); set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);
    end
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & expectTrials, csPlusTrials & rewardTrials & surpriseTrials}, 2,...
            'FluorDataField', 'ZS', 'window', [-3, 6.9], 'linespec', {'c', 'm'}); %high value, reward
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
    csPlusRewardTrials = find(csPlusTrials & rewardTrials & expectTrials);
    rn = rand(length(csPlusRewardTrials), 1);
    [~, I] = sort(rn);
    %
    subplot(2,1,1, 'FontSize', 12, 'LineWidth', 1); plot(TE.Photometry.xData, TE.Photometry.data(1).dFF(I(1:nTraces), :)', 'k'); hold on;
    phPlotAverageFromTE(TE, csPlusTrials & rewardTrials & expectTrials, 1, 'window', [-4, 7], 'linespec', {'r'});
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
    ensureFigure('all_behavior_CsPlus', 1);
    reversals = find(diff(TE.BlockNumber(csPlusTrials, :))) + 1;
    subplot(1,4,1);
    
    image(TE.Wheel.data.V(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('Velocity');
    subplot(1,4,2);
    image(TE.pupil.pupDiameterNorm(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('Pupil Diameter');    
    subplot(1,4,3); phRasterFromTE(TE, csPlusTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,4,4); phRasterFromTE(TE, csPlusTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('DAT');    
    ax = findobj(gcf, 'Type', 'Axes');
%     set(ax, 'YLim', [40 80]);
    

%% allBehavior_rasters_SFN
    firstTrial = 100;
    all_behavior_trials = csPlusTrials & rewardTrials & hitTrials & ~isnan(mean(TE.pupil.pupDiameterNorm, 2)) & TE.trialNumber >= firstTrial;
    
    savename = 'allBehavior_rasters';
    ensureFigure(savename, 1);
%     reversals = find(diff(TE.BlockNumber(csPlusTrials, :))) + 1;
    
    subplot(1,5,1); 
    eventRasterFromTE(TE, all_behavior_trials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('Licks'); ylabel('trial number');
    set(gca, 'XLim', [-4 7]);
    
%     set(gca, 'YLim', [0 max(trialCount)]);
%     set(gca, 'FontSize', 14)
    wp = [bpX2pnt(-4, 20, -4) bpX2pnt(3, 20, -4)];
    subplot(1,5,2);    
    image(TE.Wheel.data.V(all_behavior_trials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('Velocity');
    set(gca, 'YTick', []); 

    subplot(1,5,3);
    image(TE.pupil.pupDiameterNorm(all_behavior_trials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
    set(gca, 'XLim', [-4 7]);
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('Pupil Diameter');  
        set(gca, 'YTick', []); 

        
    subplot(1,5,4); phRasterFromTE(TE, all_behavior_trials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    set(gca, 'YTick', [], 'XLim', [-4 7]); 

    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,5,5); phRasterFromTE(TE, all_behavior_trials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    set(gca, 'YTick', [], 'XLim', [-4 7]); 
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('DAT');    
    ax = findobj(gcf, 'Type', 'Axes');
    formatFigurePoster([12 4], [], 20);
        if saveOn    
        saveas(gcf, fullfile(savepath, savename), 'jpeg');
        saveas(gcf, fullfile(savepath, savename), 'epsc');
        saveas(gcf, fullfile(savepath, savename), 'fig');        
        end
%     set(ax, 'YLim', [40 80]);
%% allBehavior_averages_SFN
%     ylim = [-2 8];
    savename = 'allBehavior_averages';
    h=ensureFigure(saveName, 1); 
    
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-4 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording'};
    axh = [];
    subplot(1, 5, 1); [ha, hl] = plotEventAverageFromTE(TE, csPlusTrials & rewardTrials & expectTrials, 'Port1In', varargin{:});    

% [171 55 214]/256 [237 125 49]/256
    subplot(1, 5, 4);
    [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & expectTrials}, 1,...
    'FluorDataField', 'ZS', 'window', [-4, 3], 'cmap', [171 55 214]/256, 'alpha', 0); %high value, reward
     set(gca, 'XLim', [-4, 3]);%set(gca, 'YLim', ylim);


    subplot(1, 5, 5);
    [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & expectTrials}, 2,...
        'FluorDataField', 'ZS', 'window', [-4, 3], 'cmap', [237 125 49]/256, 'alpha', 0); %high value, reward
    legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    xlabel('time from cue (s)');      set(gca, 'XLim', [-4, 3]);

    formatFigurePoster([12 2], [], 20);
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end    
    %%
    ensureFigure('all_behavior_CsMinus', 1);
    reversals = find(diff(TE.BlockNumber(csMinusTrials, :))) + 1;
    subplot(1,4,1);
    title('Velocity');
    image(TE.Wheel.data.V(csMinusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    subplot(1,4,2);
    image(TE.pupil.pupDiameterNorm(csMinusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('Pupil Diameter');    
    subplot(1,4,3); phRasterFromTE(TE, csMinusTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,4,4); phRasterFromTE(TE, csMinusTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('DAT');    
    
    
    %% plot chat, dat, and pupil averages for cue responses
    
    savename = 'glm_averages';
    ensureFigure(savename, 1); ax = axes;
    yyaxis right;  ax.YColor = [0 0 0]; % ylabel('Pupil Diameter (norm.)');
    plotPupilAverageFromTE(TE, avgTrials & ~isnan(TE.pupil_cs), 'linespec', {'k'}, 'window', [-1 3]); hold on;
    set(gca, 'YLim', [1 1.12]);
    yyaxis left; ax.YColor = [0 0 0]; hold on;
    [ha, hl] = phPlotAverageFromTE(TE, glmTrials, 1,...
    'FluorDataField', 'ZS', 'window', [-1, 3], 'cmap', [171 55 214]/256, 'alpha', 0); %high value, reward
     set(gca, 'XLim', [-4, 3]);%set(gca, 'YLim', ylim);


    [ha, hl] = phPlotAverageFromTE(TE, glmTrials, 2,...
        'FluorDataField', 'ZS', 'window', [-1, 3], 'cmap', [237 125 49]/256, 'alpha', 0); %high value, reward
    xlabel('time from cue (s)');      set(gca, 'XLim', [-1, 3]);
    formatFigurePoster([5, 3], [], 24);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    
    %% make glm (NOT really a glm though, not using any linking functions)
    chat_bl = TE.phPeakMean_baseline(1).data(glmTrials);
    pup_bl = TE.pupilBaseline(glmTrials);
    wheel_bl = TE.wheelBaseline(glmTrials);
    pup_cs = TE.pupil_cs(glmTrials); % but mouse blinks sometimes so don't include this as regressor initially
    chat_cs = TE.phPeakMean_cs(1).data(glmTrials);
    dat_cs = TE.phPeakMean_cs(2).data(glmTrials);
    table_chat = table(pup_bl, wheel_bl, chat_cs);
    
    mdl_chat = fitglm(table_chat);
    
    %% ChAT Model
    glmTrials = csPlusTrials & rewardTrials & hitTrials;
    shiftPoints = 0.3 * 20;  % 0.3s lag and 20Hz sample rate
    pupScatter = TE.pupil.pupDiameter(:,1 + shiftPoints:end);
    chatScatter = TE.Photometry.data(1).ZS(:,1:end - shiftPoints);
    ensureFigure('pup_chat_scatter', 1); scatter(pupScatter(:),chatScatter(:), '.');
    
    glmEndTime = 3;

    input_pup = TE.pupil.pupDiameterNorm(glmTrials, bpX2pnt(-1, 20, -4) + shiftPoints :bpX2pnt(glmEndTime, 20, -4) + shiftPoints);
    input_pup = input_pup'; % just concatenate the trials together
    input_pup = input_pup(:); % just concatenate the trials together
    include = ~isnan(input_pup); 

    input = TE.Photometry.data(1).ZS(glmTrials, bpX2pnt(-1, 20, -4):bpX2pnt(glmEndTime, 20, -4));    
    input = input';
    numPoints = size(input, 1);
    numTrials = size(input, 2);
    x = diag(ones(1,numPoints));
    x2 = repmat(x, numTrials, 1);
    input = input(:);
     
    input = zscore(input(include));
    input_pup = zscore(input_pup(include));
    x2 = x2(include,:);

    % regression with design matrix (time as the regressor) 
    
    B_chat = regress(input,x2);
    ensureFigure('test', 1); plot(B_chat);
    
    
    % regression with time and pupil
    
    x3 = [x2 input_pup];
    
    B1_chat = regress(input,x3);
    ensureFigure('test', 1); plot(B1_chat);
    
    %
    shuff_pup = input_pup(randperm(length(input_pup)));
    shuff_x2 = x2(randperm(size(x2,1)),randperm(size(x2,2)));
    
    nrFolds = 10;
    for iFolds = 1:nrFolds
        idx = randperm(length(input));
        idx = idx(1:round(length(input)/10));
        cIdx = false(1,length(input));
        cIdx(idx) = true;
        
        B1_chat = regress(input(~cIdx),[x2(~cIdx,:) input_pup(~cIdx)]); %full model
        B2_chat = regress(input(~cIdx),[x2(~cIdx,:) shuff_pup(~cIdx)]); %time model
        B3_chat = regress(input(~cIdx),[shuff_x2(~cIdx,:) input_pup(~cIdx)]); %pupil model
        
        Y1{iFolds} = B1_chat' * ([x2(cIdx,:) input_pup(cIdx)])';
        Y2{iFolds} = B2_chat' * ([x2(cIdx,:) shuff_pup(cIdx)])';
        Y3{iFolds} = B3_chat' * ([shuff_x2(cIdx,:) input_pup(cIdx)])';
        Y{iFolds} = input(cIdx);
        
    end
    
    fullY = cat(1,Y{:});
    fullY1 = cat(2,Y1{:})';
    fullY2 = cat(2,Y2{:})';
    fullY3 = cat(2,Y3{:})';
    
    Rsq_chat(1) = corr2(fullY,fullY1).^2;
    Rsq_chat(2) = corr2(fullY,fullY2).^2;
    Rsq_chat(3) = corr2(fullY,fullY3).^2;
    
%% DAT Model
    shiftPoints = 0.3 * 20;  % 0.3s lag and 20Hz sample rate
    pupScatter = TE.pupil.pupDiameter(:,1 + shiftPoints:end);
    datScatter = TE.Photometry.data(2).ZS(:,1:end - shiftPoints);
    ensureFigure('pup_dat_scatter', 1); scatter(pupScatter(:),datScatter(:), '.');
    

    input = TE.Photometry.data(2).ZS(glmTrials, bpX2pnt(-1, 20, -4):bpX2pnt(glmEndTime, 20, -4));
    input_pup = TE.pupil.pupDiameterNorm(glmTrials, bpX2pnt(-1, 20, -4) + shiftPoints :bpX2pnt(glmEndTime, 20, -4) + shiftPoints);
    input_pup = input_pup';
    input_pup = input_pup(:);
    include = ~isnan(input_pup);
    
    input = input';
    numPoints = size(input, 1);
    numTrials = size(input, 2);
    x = diag(ones(1,numPoints));
    x2 = repmat(x, numTrials, 1);
    input = input(:);
     
    input = zscore(input(include));
    input_pup = zscore(input_pup(include));
    x2 = x2(include,:);

    % regression with design matrix (time as the regressor) 
    
    B_dat = regress(input,x2);
    ensureFigure('test', 1); plot(B_dat);
    
    
    % regression with time and pupil
    
    x3 = [x2 input_pup];
    
    B1_dat = regress(input,x3);
    ensureFigure('test', 1); plot(B1_dat);
    
    %
    shuff_pup = input_pup(randperm(length(input_pup)));
    shuff_x2 = x2(randperm(size(x2,1)),randperm(size(x2,2)));
    
    nrFolds = 10;
    for iFolds = 1:nrFolds
        idx = randperm(length(input));
        idx = idx(1:round(length(input)/10));
        cIdx = false(1,length(input));
        cIdx(idx) = true;
        
        B1_dat = regress(input(~cIdx),[x2(~cIdx,:) input_pup(~cIdx)]); %full model
        B2_dat = regress(input(~cIdx),[x2(~cIdx,:) shuff_pup(~cIdx)]); %time model
        B3_dat = regress(input(~cIdx),[shuff_x2(~cIdx,:) input_pup(~cIdx)]); %pupil model
        
        Y1{iFolds} = B1_dat' * ([x2(cIdx,:) input_pup(cIdx)])';
        Y2{iFolds} = B2_dat' * ([x2(cIdx,:) shuff_pup(cIdx)])';
        Y3{iFolds} = B3_dat' * ([shuff_x2(cIdx,:) input_pup(cIdx)])';
        Y{iFolds} = input(cIdx);
        
    end
    
    fullY = cat(1,Y{:});
    fullY1 = cat(2,Y1{:})';
    fullY2 = cat(2,Y2{:})';
    fullY3 = cat(2,Y3{:})';
    
    Rsq_dat(1) = corr2(fullY,fullY1).^2;
    Rsq_dat(2) = corr2(fullY,fullY2).^2;
    Rsq_dat(3) = corr2(fullY,fullY3).^2;    
    
    %% plot beta weights for ChAT and DAT
    
    
    xdata = [linspace(-1, 3, 81) 4];
    savename = 'glm_betaWeights_both';
    ensureFigure(savename, 1);
    subplot(1,2,1);
    hl = plot(xdata, B1_chat); hold on
    hl.Color = [171 55 214]/256; hl.LineWidth = 2;
        set(gca, 'XLim', [-1 4.5]);
%     formatFigurePoster([5, 3], [], 24);
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end

%     savename = 'glm_betaWeights_dat';
%     ensureFigure(savename, 1);
subplot(1,2,2);
hl = plot(xdata, B1_dat); hold on
    hl.Color = [237 125 49]/256; hl.LineWidth = 2;
    set(gca, 'XLim', [-1 4.5]);
    formatFigurePoster([9, 3], [], 24);
    
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    

%% bar graphs of R2 values

    savename = 'glm_Rsq_both';
    ensureFigure(savename, 1);
    subplot(1,2,1);
    bar(Rsq_chat, 'FaceColor', [171 55 214]/256);
    set(gca, 'YLim', [0 0.4], 'XTickLabel', {'full', 'time', 'pupil'});
    ylabel('R squared');
    subplot(1,2,2);
    bar(Rsq_dat, 'FaceColor', [237 125 49]/256);    
    set(gca, 'YLim', [0 0.4], 'XTickLabel', {'full', 'time', 'pupil'}, 'YTick', []);
    
    formatFigurePoster([9, 3], [], 24);
    
if saveOn
    saveas(gcf, fullfile(savepath, [savename '.fig']));
    saveas(gcf, fullfile(savepath, [savename '.jpg']));   
    saveas(gcf, fullfile(savepath, [savename '.epsc']));   
end    
    