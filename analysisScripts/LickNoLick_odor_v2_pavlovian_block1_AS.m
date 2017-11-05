
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
    BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
    BL{end + 1} = [2 4];    
end


    

TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL);
%%


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
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = uigetdir;
sep = strfind(TE.filename{1}, '.');
subjectName = TE.filename{1}(1:sep(1)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

%% grouped sessions
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
        try
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
end
%% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks    blocks 2 and 3
    validTrials = filterTE(TE, 'reject', 0);
    Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
    rewardTrials = filterTE(TE, 'trialType', [1 3]);
    
    cuedReward = filterTE(TE, 'trialType', 1, 'reject', 0);
    omit = filterTE(TE, 'trialType', 2, 'reject', 0);
    uncuedReward = filterTE(TE, 'trialType', 3, 'reject', 0);

    trialCount = [1:length(TE.filename)]';
    
%%
%% photometry averages, zscored
    saveName = 'phAvgs';  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);

    pm = [1 2];
    
    % - 6 0 4
    if ismember(1, channels)
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {cuedReward, omit, uncuedReward}, 1,...
            'FluorDataField', 'ZS', 'window', [-4, 5.5], 'linespec', {'b', 'k', 'c'}); %high value, reward
        legend(hl, {'cued', 'omit', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);
    end
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {cuedReward, omit, uncuedReward}, 2,...
            'FluorDataField', 'ZS', 'window', [-4, 5.5], 'linespec', {'b', 'k', 'c'}); %high value, reward
        legend(hl, {'cued', 'omission', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        ylabel('VTA dF/F Zscored'); xlabel('time from cue (s)'); 
    end
    
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end    
    
    
%% Raster Plots
CLimFactor = 2;


for channel = channels

    saveName = ['rasters_ch' num2str(channel)];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);
    

%     prcd = TE.Photometry.data(1).dFF(prt, :);
    subplot(1,4,1); 
    eventRasterFromTE(TE, cuedReward, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('cued'); ylabel('trial number');
    set(gca, 'YLim', [0 sum(cuedReward)]);
    set(gca, 'XLim', [-4 7]); 
    set(gca, 'FontSize', 14)
   

    
    subplot(1,4,2); phRasterFromTE(TE, cuedReward, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 14);

    subplot(1,4,3); phRasterFromTE(TE, omit, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 14); title('omission');

    subplot(1,4,4); phRasterFromTE(TE, uncuedReward, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        set(gca, 'FontSize', 14); title('uncued');       

    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
end

%%

pm = [4 4];
saveName = 'examples';
ensureFigure(saveName, 1);
for counter = 1:16
    subplot(pm(1), pm(2), counter);
    plot(TE.Photometry.xData, TE.Photometry.data(1).dFF(counter, :));
end

    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
    