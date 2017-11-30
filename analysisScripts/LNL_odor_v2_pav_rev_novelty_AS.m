% 4/10/17  Analysis script for novel odorant learning combined with pavlovian reversals using LickNoLick_Odor_V2
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
% basepath = uigetdir;
basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_Novel_Odorants\';
% basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_20and17_combined\';
% if length(unique(TE.filename)) > 1
%     sep = strfind(TE.filename{1}, '_');
%     subjectName = TE.filename{1}(1:sep(2)-1);
% else
    sep = strfind(TE.filename{1}, '.');
    subjectName = TE.filename{1}(1:sep(1)-1);
% end
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);


%%
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'ReinforcementOutcome', 'Reward'));
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
% NOTE!  updated these lookups to be general across blocks (don't make use
% of trial type, etc.)
    validTrials = filterTE(TE, 'reject', 0);
    Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
    Odor2Trials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);    
    Odor3Trials = filterTE(TE, 'OdorValveIndex', 3, 'reject', 0);  
    
    rewardTrials = filterTE(TE, 'ReinforcementOutcome', 'Reward', 'reject', 0);
    hitTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    missTrials = filterTE(TE, 'trialOutcome', -1, 'reject', 0);
    punishTrials = filterTE(TE, 'ReinforcementOutcome', 'Punish', 'reject', 0);    
    neutralTrials = filterTE(TE, 'ReinforcementOutcome', 'Neutral', 'reject', 0);
    csPlusTrials = filterTE(TE, 'CSValence', 1, 'reject', 0);
    csMinusTrials = filterTE(TE, 'CSValence', -1, 'reject', 0);
    uncuedReward = filterTE(TE, 'ReinforcementOutcome', 'Reward', 'OdorValveIndex', 0, 'reject', 0);
    uncuedPunish = filterTE(TE, 'ReinforcementOutcome', 'Punish', 'OdorValveIndex', 0,'reject', 0);


    trialCount = [1:length(TE.filename)]';
%% 


%% PH Rasters by odor,
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
    subplot(1,6,1); 
    eventRasterFromTE(TE, Odor1Trials, 'Port1In', 'trialNumbering', 'global',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('Odor 1'); ylabel('trial number');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    

    
%     subplot(1,6,2); phRasterFromTE(TE, csPlusTrials, channel, 'CLimFactor', CLimFactor, 'CLim', CLim{channel}, 'trialNumbering', 'global');
%         set(gca, 'FontSize', 14)
    subplot(1,6,2); phRasterFromTE(TE, Odor1Trials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
        set(gca, 'FontSize', 14)
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines        
        
    subplot(1,6,3); 
    eventRasterFromTE(TE, Odor2Trials, 'Port1In', 'trialNumbering', 'global',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    title('Odor 2');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', [0 max(trialCount)]);
    xlabel('time from cue (s)');     
    set(gca, 'FontSize', 14)
    
    subplot(1,6,4); phRasterFromTE(TE, Odor2Trials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    set(gca, 'FontSize', 14)
    
    subplot(1,6,5); 
    eventRasterFromTE(TE, Odor3Trials, 'Port1In', 'trialNumbering', 'global',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    title('Odor 3');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    
    subplot(1,6,6); phRasterFromTE(TE, Odor3Trials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    set(gca, 'FontSize', 14)    
    colormap('parula');  

    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
end

%%
%% PH Rasters novel odor focus, both channels on same plot
% right now I'm only considering case where I've introduced a novel odor
% once
firstNTrials = 20;
CLimFactor = 2;
trialStart = TE.Photometry.xData(1);
reversals = find(TE.BlockChange);
trialsToShow = find(Odor3Trials);
trialsToShow = trialsToShow(1:min(firstNTrials, length(trialsToShow)));
    saveName = [subjectName '_novelOdor'];
    h=ensureFigure(saveName, 1);
    mcLandscapeFigSetup(h);
    

%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,3,1); 
    eventRasterFromTE(TE, trialsToShow, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('Novel odor (odor index = 3)'); ylabel('trial number');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', [0 length(trialsToShow)]);
    set(gca, 'FontSize', 14)    

    
    subplot(1,3,2); phRasterFromTE(TE, trialsToShow, 1, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 14); title('ChAT');
    
    subplot(1,3,3); phRasterFromTE(TE, trialsToShow, 2, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 14); title('DAT');
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end        

%% compare learning of reversed and novel odorants for ChAT and DAT 
saveName = [subjectName '_longitudinal_novelOdorant'];
h=ensureFigure(saveName, 1); 
mcLandscapeFigSetup(h);
newCsPlusTrials = csPlusTrials & (Odor1Trials | Odor2Trials) & (TE.BlockNumber == 8 | TE.BlockNumber == 9);
subplot(3,1,1); plot(TE.csLicks.rate(Odor3Trials),'go', 'LineStyle', 'none');  hold on
plot(TE.csLicks.rate(newCsPlusTrials),'r+', 'LineStyle', 'none');
% subplot(4,1,1); scatter(1:length(allTE.epoch), cuedReward.anticipatoryLicks2);
% maxLR = max(TE.csLicks.rate(csPlusTrials));
% % subplot(4,1,1); plot(smooth(cuedReward.anticipatoryLicks1 / mean(cuedReward.blLicks)));
% % maxLR = smooth(max(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
% hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxLR, 'g', 'Marker', 'none');
% hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxLR, 'r', 'Marker', 'none');
ylabel('antic. licks'); title('CS+ and/or reward trials, dF/F is Zscored');
% set(gca, 'XLim', [1 length(trialCount)]);

if ismember(1, channels)
    subplot(3,1,2); plot(TE.phPeakPercentile_cs(1).data(Odor3Trials), 'g.', 'LineStyle', 'none'); hold on;
    plot(TE.phPeakPercentile_cs(1).data(newCsPlusTrials), 'r.', 'LineStyle', 'none');
    ylabel('BF: Cue ZS');
%     set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(2, channels)
    subplot(3,1,3); plot(TE.phPeakPercentile_cs(2).data(Odor3Trials), 'g.', 'LineStyle', 'none'); hold on;
    plot(TE.phPeakPercentile_cs(2).data(newCsPlusTrials), 'r.', 'LineStyle', 'none');
    ylabel('VTA: Cue ZS');
%     set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end
    