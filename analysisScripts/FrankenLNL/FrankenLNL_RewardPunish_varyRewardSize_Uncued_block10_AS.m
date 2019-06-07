%% for rewardPunishBlocks == 10 (uncued)
% trial type 1 = 10uL, trial type 2 = 5uL, trial type 3 = 2uL, trial type 4
% = 0uL

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

duration = length(TE.Photometry.xData) / TE.Photometry.sampleRate;

%% basepath = uigetdir;
basepath = 'Z:\SummaryAnalyses\Franken_LNL_rewardSize';
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

%%
TE.licks_us = countEventFromTE(TE, 'Port1In', usWindow, TE.Us);
usWindow = [0 .5];  
minRewardLickRate = percentile(TE.licks_us.rate(rewardTrials), 0.2);
frankenLNL_conditions;




bigHit = trialsByType{1} & (TE.licks_us.rate >= minRewardLickRate);
mediumHit = trialsByType{2} & (TE.licks_us.rate >= minRewardLickRate);
smallHit = trialsByType{3} & (TE.licks_us.rate >=minRewardLickRate);

bigMiss = trialsByType{1} & (TE.licks_us.rate < minRewardLickRate);
mediumMiss = trialsByType{2} & (TE.licks_us.rate < minRewardLickRate);
smallMiss = trialsByType{3} & (TE.licks_us.rate < minRewardLickRate);



%%
saveName = 'Averages';
ensureFigure(saveName, 1);

trials = {trialsByType{1}, trialsByType{2}, trialsByType{3}, neutralTrials};
window = [-0.5 3];
PhotometryField = 'Photometry';
subplot(1,3,1);
        [ha, hl] = phPlotAverageFromTE(TE, trials, 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'b', 'c', 'm', 'k'}, 'PhotometryField', PhotometryField); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [-0.1 0.1]);
        set(gca, 'XLim', window); ylabel('Fluor ZS');
        title('Ch 1');
        
subplot(1,3,2);
        [ha, hl] = phPlotAverageFromTE(TE, trials, 2, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'b', 'c', 'm', 'k'}, 'PhotometryField', PhotometryField); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [-0.1 0.1]);
        set(gca, 'XLim', window); ylabel('Fluor ZS');        
        title('Ch 2');
subplot(1,3,3);
        varargin = {'window', window, 'zeroTimes', TE.Us, 'window', window, 'linespec', {'b', 'c', 'm', 'k'}};
        axh = [];
        [ha, hl] = plotEventAverageFromTE(TE, trials, 'Port1In', varargin{:});  
        addStimulusPatch(gca, [-0.1 0.1]);
        ylabel('licks (1/s)'); xlabel('time from cue(s)');  set(gca, 'XLim', window);
        legend(hl, {'10', '5', '2', 'none'}, 'Box', 'off', 'Location', 'best'); 
        title('Licks');
        
        
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end        
        
%% hit vs miss

saveName = 'Averages_hit_vs_miss';
ensureFigure(saveName, 1);


window = [-0.5 3];
PhotometryField = 'PhotometryExpFit';
subplot(3,2,1);
        [ha, hl] = phPlotAverageFromTE(TE, {bigHit, bigMiss}, 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'g', 'r'}, 'PhotometryField', PhotometryField); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'best'); legend('boxoff');
        addStimulusPatch(gca, [-0.1 0.1]);
        set(gca, 'XLim', window); ylabel('Big ZS');
        title('Ch 1');
        
subplot(3,2,2);
        [ha, hl] = phPlotAverageFromTE(TE, {bigHit, bigMiss}, 2, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'g', 'r'}, 'PhotometryField', PhotometryField); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [-0.1 0.1]);
        set(gca, 'XLim', window); 
        title('Ch 2');        
        
subplot(3,2,3);
        [ha, hl] = phPlotAverageFromTE(TE, {mediumHit, mediumMiss}, 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'g', 'r'}, 'PhotometryField', PhotometryField); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [-0.1 0.1]);
        set(gca, 'XLim', window); ylabel('Medium ZS');        
        
subplot(3,2,4);
        [ha, hl] = phPlotAverageFromTE(TE, {mediumHit, mediumMiss}, 2, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'g', 'r'}, 'PhotometryField', PhotometryField); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [-0.1 0.1]);
        set(gca, 'XLim', window); 
        
        
        
        
subplot(3,2,5);
        [ha, hl] = phPlotAverageFromTE(TE, {smallHit, smallMiss}, 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'g', 'r'}, 'PhotometryField', PhotometryField); %high value, reward
        addStimulusPatch(gca, [-0.1 0.1]);
        set(gca, 'XLim', window); ylabel('Small ZS');        
        
subplot(3,2,6);
        [ha, hl] = phPlotAverageFromTE(TE, {smallHit, smallMiss}, 2, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'g', 'r'}, 'PhotometryField', PhotometryField); %high value, reward
        addStimulusPatch(gca, [-0.1 0.1]);
        set(gca, 'XLim', window); 
        
        
        
