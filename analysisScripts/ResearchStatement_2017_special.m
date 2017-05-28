saveOn = 0;
%% 
sessions = bpLoadSessions;
% ***Loaded DC_17_lickNoLick_Odor_v2_Apr09_2017_Session1.mat from Z:\FitzRig2\Data\DC_17\lickNoLick_Odor_v2\Session Data\***
%%
TE = makeTE_LNL_odor_V2(sessions);
%%
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRGECO1a (ch2)
channels=[]; dFFMode = {}; BL = {};
if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'expFit';
    BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
    BL{end + 1} = [2 4];    
end


    

TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL, 'downsample', 305);
%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = uigetdir;
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);
%% extract peak trial dFF responses to cues and reinforcement and lick counts
% csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);
csWindow = zeros(nTrials, 2);
csWindow(:,2) = cellfun(@(x,y,z) max(x(end), y(end)) - z(1), TE.AnswerLick, TE.AnswerNoLick, TE.Cue); 
% max 1) to select AnswerNoLick time stamp for no lick trials (unused state contains NaN)
% 2) to select AnswerLick time stamp for lick trials (AnswerLick follows
% AnswerNoLick state)

ch1CsWindow = [0.5 1.5];
ch2CsWindow = [0.5 1.5];

TE.csLicks = countEventFromTE(TE, 'Port1In', csWindow, TE.Cue);

usWindow = [0 0.5];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'


for channel = channels
    TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'mean', 'phField', 'ZS');
    if channel == 1
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.9, 'phField', 'ZS');
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');
    elseif channel == 2
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.9, 'phField', 'ZS');        
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');  
    end
end



TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);


%%
truncateSessionsFromTE(TE, 'init');

%% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks    blocks 2 and 3
    validTrials = filterTE(TE, 'reject', 0);
%     Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
%     Odor2Trials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);    
%     
%     rewardTrials = filterTE(TE, 'trialType', 1, 'reject', 0);
%     hitTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
%     missTrials = filterTE(TE, 'trialOutcome', -1, 'reject', 0);
%     punishTrials = filterTE(TE, 'trialType', 3, 'reject', 0);    
%     neutralTrials = filterTE(TE, 'trialType', [2 4], 'reject', 0);
%     block2Trials = filterTE(TE, 'BlockNumber', 2, 'reject', 0);
%     block3Trials = filterTE(TE, 'BlockNumber', 3, 'reject', 0);
    csPlusTrials = filterTE(TE, 'trialType', [1 2], 'reject', 0);
    csMinusTrials = filterTE(TE, 'trialType', [3 4], 'reject', 0);    
    trialTypes = 1:4;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end




%% PH Rasters, CS+, CS-
CLimFactor = 2;
CLimFactor2 = 2.5;
trialRange = [20 80];
reversals = find(diff(TE.BlockNumber(csPlusTrials)));

    saveName = ['researchStatement_reversals_phRasters_dualChannel'];
    h=ensureFigure(saveName, 1);
%     mcPortraitFigSetup(h);
    

%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,3,1); 
    eventRasterFromTE(TE, csPlusTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
%     title('CS+'); ylabel('trial number');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', trialRange, 'TickDir', 'Out');
    set(gca, 'FontSize', 10)
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    set(gca, 'XLim', [-2 6]);
    
    subplot(1,3,2); phRasterFromTE(TE, csPlusTrials, 1, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 10, 'TickDir', 'Out')
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines        
            set(gca, 'YLim', trialRange);set(gca, 'XLim', [-2 6]);
     
    subplot(1,3,3); phRasterFromTE(TE, csPlusTrials, 2, 'CLimFactor', CLimFactor2, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 10, 'TickDir', 'Out')
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines       
        set(gca, 'YLim', trialRange);set(gca, 'XLim', [-2 6]);
  


    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
    end

%%

%% hit vs miss averages, zscored
    % - 6 0 4
    saveName = [subjectName '_hitVSmiss_ch1'];  
    h=ensureFigure(saveName, 1); 
    axes('FontSize', 10, 'LineWidth', 1); 
    [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 1,...
    'FluorDataField', 'ZS', 'window', [-2, 6], 'linespec', {'g', 'r'}); %high value, reward
    legend(hl, {'hit', 'miss'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
     set(gca, 'XLim', [-2, 6]);
    formatFigureGRC;
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));           
    end
    
    saveName = [subjectName '_hitVSmiss_ch2'];
    h=ensureFigure(saveName, 1); 
    subplot(1, 1, 1, 'FontSize', 10, 'LineWidth', 1); 
    [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 2,...
        'FluorDataField', 'ZS', 'window', [-2, 6], 'linespec', {'g', 'r'}); %high value, reward
    legend(hl, {'hit', 'miss'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    xlabel('time from cue (s)');     set(gca, 'XLim', [-2, 6]);
    formatFigureGRC;
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));           
    end    
    
    
    %% hit vs miss, licking
    
        % - 6 0 4
    saveName = [subjectName '_hitVSmiss_licking'];  
    h=ensureFigure(saveName, 1); 
    axes('FontSize', 10, 'LineWidth', 1); 
    
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-2 6], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'linespec', {'g', 'r'}};
    axh = [];
    subplot(1, 1, 1); [ha, hl] = plotEventAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 'Port1In', varargin{:});
    
    legend(hl, {'hit', 'miss'}, 'Location', 'northwest', 'FontSize', 10); legend('boxoff');
     set(gca, 'XLim', [-2, 6]);
    formatFigureGRC;
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));           
    end
    
    