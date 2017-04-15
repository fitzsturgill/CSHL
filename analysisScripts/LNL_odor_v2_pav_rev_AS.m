% 4/10/17  Analysis script for pavlovian reversals using LickNoLick_Odor_V2
% protocol



saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_LNL_odor_V2(sessions);
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRGECO1a (ch2)
channels=[]; dFFMode = {};
if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'expFit';
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
end


    

TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels);


%% extract peak trial dFF responses to cues and reinforcement and lick counts
% csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);
csWindow = zeros(nTrials, 2);
csWindow(:,2) = cellfun(@(x,y,z) max(x(end), y(end)) - z(1), TE.AnswerLick, TE.AnswerNoLick, TE.Cue); 
% max 1) to select AnswerNoLick time stamp for no lick trials (unused state contains NaN)
% 2) to select AnswerLick time stamp for lick trials (AnswerLick follows
% AnswerNoLick state)



TE.csLicks = countEventFromTE(TE, 'Port1In', csWindow, TE.Cue);

usWindow = [0 0.5];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'

for channel = channels
    TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0.5 1.5], TE.Cue, 'method', 'mean');
    TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0.5 1.5], TE.Cue, 'method', 'percentile', 'percentile', 0.8);
    TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'mean');
    if channel == 1
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'percentile', 'percentile', 0.8);
    elseif channel == 2
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'percentile', 'percentile', 0.5);  
    end
end



TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);

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


%%
truncateSessionsFromTE(TE, 'init');
%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% cross sessions bleaching curve and dual exponential fits
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
            plot(TE.Photometry.bleachFit(counter, channel).trialTemplate, 'k'); hold on;
            plot(TE.Photometry.bleachFit(counter, channel).trialFit, 'r');
        %     title(num2str(counter));    
        end
        if saveOn
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.fig'));
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.jpg'));
        end
    end
end
%% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks    blocks 2 and 3
    validTrials = filterTE(TE, 'reject', 0);
    Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
    Odor2Trials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);    
    
    rewardTrials = filterTE(TE, 'trialType', 1, 'reject', 0);
    hitTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    missTrials = filterTE(TE, 'trialOutcome', -1, 'reject', 0);
    punishTrials = filterTE(TE, 'trialType', 3, 'reject', 0);    
    neutralTrials = filterTE(TE, 'trialType', [2 4], 'reject', 0);
    block2Trials = filterTE(TE, 'BlockNumber', 2, 'reject', 0);
    block3Trials = filterTE(TE, 'BlockNumber', 3, 'reject', 0);
    csPlusTrials = filterTE(TE, 'trialType', [1 2], 'reject', 0);
    csMinusTrials = filterTE(TE, 'trialType', [3 4], 'reject', 0);    
    trialTypes = 1:4;
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
% subplot(4,1,1); plot(smooth(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
% maxLR = smooth(max(cuedReward.anticipatoryLicks2 / mean(cuedReward.blLicks)));
hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxLR, 'g', 'Marker', 'none');
hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxLR, 'r', 'Marker', 'none');
ylabel('antic. licks'); title('CS+ and/or reward trials');
set(gca, 'XLim', [1 length(trialCount)]);

if ismember(1, channels)
    subplot(5,1,2); scatter(trialCount(csPlusTrials), TE.phPeakPercentile_cs(1).data(csPlusTrials), '.');
    maxP = max(TE.phPeakPercentile_cs(1).data(csPlusTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
    ylabel('BF: CS DF/F (80%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(2, channels)
    subplot(5,1,3); scatter(trialCount(csPlusTrials), TE.phPeakPercentile_cs(2).data(csPlusTrials), '.');
    maxP = max(TE.phPeakPercentile_cs(2).data(csPlusTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
ylabel('VTA:CS DF/F (80%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(1, channels)
    subplot(5,1,4); scatter(trialCount(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), '.');
    maxP = max(TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
    ylabel('BF: US DF/F (80%)');
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
ylabel('antic. licks'); title('CS- and/or punish trials');
set(gca, 'XLim', [1 length(trialCount)]);

if ismember(1, channels)
    subplot(5,1,2); scatter(trialCount(csMinusTrials), TE.phPeakPercentile_cs(1).data(csMinusTrials), '.');
    maxP = max(TE.phPeakPercentile_cs(1).data(csMinusTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
    ylabel('BF: CS DF/F (80%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(2, channels)
    subplot(5,1,3); scatter(trialCount(csMinusTrials), TE.phPeakPercentile_cs(2).data(csMinusTrials), '.');
    maxP = max(TE.phPeakPercentile_cs(2).data(csMinusTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
ylabel('VTA:CS DF/F (80%)');
    set(gca, 'XLim', [1 length(trialCount)]); %set(gca, 'YLim', [-0.2 0.2]);
end

if ismember(1, channels)
    subplot(5,1,4); scatter(trialCount(csMinusTrials & punishTrials), TE.phPeakPercentile_us(1).data(csMinusTrials & punishTrials), '.');
    maxP = max(TE.phPeakPercentile_us(1).data(csMinusTrials & punishTrials));
    hold on; stem(trialCount(TE.BlockChange ~= 0), TE.BlockChange(TE.BlockChange ~= 0) * maxP, 'g', 'Marker', 'none');
    hold on; stem(trialCount(TE.sessionChange ~= 0), TE.sessionChange(TE.sessionChange ~= 0) * maxP, 'r', 'Marker', 'none');
    ylabel('BF: US DF/F (80%)');
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

    
    subplot(1,4,2); phRasterFromTE(TE, csPlusTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
        set(gca, 'FontSize', 14)
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines        
        
    subplot(1,4,3); 
    eventRasterFromTE(TE, csMinusTrials, 'Port1In', 'trialNumbering', 'global',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    title('CS-');
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)
    
    subplot(1,4,4); phRasterFromTE(TE, csMinusTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
    set(gca, 'FontSize', 14)
    xlabel('time from cue (s)');         


    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
end

%% photometry averages
    saveName = [subjectName '_phAvgs'];  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);

    pm = [2 2];
    
    % - 6 0 4
    if ismember(1, channels)
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, neutralTrials}, 1, 'window', [3, 7], 'linespec', {'b', 'r', 'k'}); %high value, reward
        legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F'); textBox(subjectName);
    end
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, neutralTrials}, 2, 'window', [3, 7], 'linespec', {'b', 'r', 'k'}); %high value, reward
        legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        ylabel('VTA dF/F'); xlabel('time from cue (s)'); 
    end
    
    % - 6 0 4
    if ismember(1, channels)    
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 1, 'window', [-4, 7], 'linespec', {'c', 'm'}); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('CS+, outcomes'); set(gca, 'XLim', [-4, 7]);
    end
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 2, 'window', [-4, 7], 'linespec', {'c', 'm'}); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        xlabel('time from cue (s)');     set(gca, 'XLim', [-4, 7]);
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
        title('CS+, Cue'); ylabel('BF dF/F'); xlabel('VTA dF/F'); textBox(subjectName);
    end

    if all(ismember([1 2], channels))    
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1);
        scatter(TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), '.'); hold on;
        fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
        fob = fit(TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), 'poly1', fo); 
        plot(fob,'predfunc'); legend off;    
        title('CS+, Reward'); ylabel('BF dF/F'); xlabel('VTA dF/F');  
    end
    
    if ismember(1, channels)    
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1);
        scatter(TE.csLicks.rate(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), 'o'); hold on;
        fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
        try
            fob = fit(TE.csLicks.rate(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(1).data(csPlusTrials & rewardTrials), 'poly1', fo); 
            plot(fob,'predfunc'); legend off;      
        end
        title('CS+, reward vs. cue licks'); ylabel('BF dF/F'); xlabel('Antic. Lick Rate');    
    end
    
    if ismember(2, channels)        
        subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1);
        scatter(TE.csLicks.rate(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials), 'o'); hold on;
        fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
        try
            fob = fit(TE.csLicks.rate(csPlusTrials & rewardTrials), TE.phPeakPercentile_us(2).data(csPlusTrials & rewardTrials), 'poly1', fo); 
            plot(fob,'predfunc'); legend off;         
        end
        ylabel('VTA dF/F'); xlabel('Antic. Lick Rate');
    end
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end    

 %% single trial traces
    nTraces = 50;
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

