%% shujing's analysis script

% what is going on
saveOn = 0;
%%
saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_LNL_odor_V2(sessions);

%%
% assume that photometry channels are consistent across sessions
channels=[]; dFFMode = {}; BL = {};
if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
% 'expFit' subtracts biexponential fit to initial bleaching transient within
% trials (flattens this artifact), also try 'simple'
    dFFMode{end+1} = 'expFit'; 
    BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
    BL{end + 1} = [1 4];    
end
  
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL);

%% extract peak trial dFF responses to cues and reinforcement and lick counts
% zero is defined as time of cue- see call to
% processTrialAnalysis_Photometry2
nTrials = length(TE.filename);

csWindow = [0.5 1.5];
TE.csLicks = countEventFromTE(TE, 'Port1In', csWindow, TE.Cue);

usWindow = [0 0.75];
% line below extracts time when the us occurs (accounting for the various
% outcome state possibilities)
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros); %wider window for counting US licks than photometry US response

for channel = channels
    TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, BL{channel}, [], 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
    TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.9, 'phField', 'ZS');
    TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'percentile', 'percentile', 0.9, 'phField', 'ZS');
end

%% save data in a base directory, code below creates a folder named according to subject (e.g. DAT_1) and sets the save path within
basepath = uigetdir;
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);


%% exclude trials at end of session where the mouse stops licking
rewardTrialsTrunc = filterTE(TE, 'trialType', [1 3]);
truncateSessionsFromTE(TE, 'init', 'usLicks', rewardTrialsTrunc);
% left/right arrow to adjust truncation point.  up/down arrow to switch
% sessions 'u' to update


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
            saveas(gcf, fullfile(savepath, [figname '.fig']));
            saveas(gcf, fullfile(savepath, [figname '.jpg']));
        end
    catch
    end
end
%% generate trial lookups for different combinations of conditions
% see Pavlovian_reversals_blocks    blocks 2
    validTrials = filterTE(TE, 'reject', 0);
    Odor1Trials = filterTE(TE, 'OdorValveIndex', 1, 'reject', 0);
    Odor2Trials = filterTE(TE, 'OdorValveIndex', 2, 'reject', 0);    
    
    rewardTrials = filterTE(TE, 'trialType', [1], 'reject', 0);
    punishTrials = filterTE(TE, 'trialType', [3], 'reject', 0);
    
    cuedReward = filterTE(TE, 'OdorValveIndex', 1, 'ReinforcementOutcome', 'Reward', 'reject', 0);
    uncuedReward = filterTE(TE, 'OdorValveIndex', 0, 'ReinforcementOutcome', 'Reward', 'reject', 0);
    rewardOmission = filterTE(TE, 'OdorValveIndex', 1, 'ReinforcementOutcome', 'Neutral', 'reject', 0);
    cuedPunish = filterTE(TE, 'OdorValveIndex', 2, 'ReinforcementOutcome', 'Punish', 'reject', 0);
    punishOmission = filterTE(TE, 'OdorValveIndex', 2, 'ReinforcementOutcome', 'Neutral', 'reject', 0);
    
    % NOTE INCONSISTENCY ACROSS BLOCK 1 AND 2 OF TRIAL TYPES-  IN BLOCK 2
    % TRIAL TYPE 3 = CUED PUNISHMENT WHEREAS IN BLOCK 1 = UNCUED REWARD
    trialTypes = 1:4; % max 6 trial types when using shujing_blocks, block #2

    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end

    trialCount = [1:length(TE.filename)]';
    
    % example of how to select subset of trials- select trials from
%     range 1-100
    % 
%     cuedReward_subset = logical(zeros(size(cuedReward)));
%     cuedReward_subset(1:100) = cuedReward(1:100);
    
    
    %% lick and photometry rasters
    
    
    CLimFactor = 2; % scale each image +/-  2 * global baseline standard deviation

for channel = channels

    saveName = [subjectName '_rasters_ch' num2str(channel)];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);
    
    subplot(1,3,1); % lick raster
    eventRasterFromTE(TE, cuedReward, 'Port1In', 'trialNumbering', 'global',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('Cued Reward'); ylabel('trial number');
    set(gca, 'XLim', [-4 6]); 
%     set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)        
    
    % photometry raster
    subplot(1,3,2); phRasterFromTE(TE, cuedReward, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
        set(gca, 'FontSize', 14)

    subplot(1,3,3); phRasterFromTE(TE, cuedPunish, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
    title('Cued Punish');
    set(gca, 'FontSize', 14)  
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
%         saveas(gcf, fullfile(savepath, [saveName '.epsc']));   % for illustrator, uncomment        
    end

end

%% Averages, Reward
    saveName = [subjectName '_rewardAvgs'];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);

    pm = [3 1]; % subplot matrix
    
    % - 6 0 4
    
        % lick averages
    % I don't even remember why I supply trialNumbering parameter right
    % now...
    varargin = {'trialNumbering', 'consecutive',...
        'binWidth', 0.5, 'window', [-4 6], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'linespec', {'b', 'c', 'k'}};
    axh = [];
    subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, {cuedReward, uncuedReward, rewardOmission}, 'Port1In', varargin{:});
    legend(hl, {'cued', 'uncued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Licks'); ylabel('licks (s)'); xlabel('time from cue (s)'); 
    
    
    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
    % channel provided as the 3rd argument
    [ha, hl] = phPlotAverageFromTE(TE, {cuedReward, uncuedReward, rewardOmission}, 1,...
        'FluorDataField', 'dFF', 'window', [-4, 6], 'linespec', {'b', 'c', 'k'}); %high value, reward
    legend(hl, {'cued', 'uncued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Ch1'); ylabel('Ch1 dF/F');
                        
    if ismember(2, channels)
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        % channel provided as the 3rd argument
        [ha, hl] = phPlotAverageFromTE(TE, {cuedReward, uncuedReward, rewardOmission}, 2,...
            'FluorDataField', 'dFF', 'window', [-4, 6], 'linespec', {'b', 'c', 'k'}); %high value, reward
        legend(hl, {'cued', 'uncued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Ch2'); ylabel('Ch2 dF/F');               
    end
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));   % for illustrator, uncomment        
    end    

%% Averages, Punish
    saveName = [subjectName '_punishAvgs'];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);

    pm = [3 1]; % subplot matrix
    
    % - 6 0 4
    
        % lick averages
    % I don't even remember why I supply trialNumbering parameter right
    % now...
    varargin = {'trialNumbering', 'consecutive',...
        'binWidth', 0.5, 'window', [-4 6], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'linespec', {'r', 'k'}};
    axh = [];
    subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, {cuedPunish, punishOmission}, 'Port1In', varargin{:});
    legend(hl, {'cued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Licks'); ylabel('licks (s)'); xlabel('time from cue (s)'); 
    
    
    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
    % channel provided as the 3rd argument
    [ha, hl] = phPlotAverageFromTE(TE, {cuedPunish, punishOmission}, 1,...
        'FluorDataField', 'dFF', 'window', [-4, 6], 'linespec', {'r', 'k'}); %high value, reward
    legend(hl, {'cued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Ch1'); ylabel('Ch1 dF/F');
                        
    if ismember(2, channels)
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        % channel provided as the 3rd argument
        [ha, hl] = phPlotAverageFromTE(TE, {cuedPunish, punishOmission}, 2,...
            'FluorDataField', 'dFF', 'window', [-4, 6], 'linespec', {'r', 'k'}); %high value, reward
        legend(hl, {'cued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Ch2'); ylabel('Ch2 dF/F');             
    end
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));   % for illustrator, uncomment        
    end        

    
    %%
        %% summary statistics
    cComplete_summary = struct(...
        'cuedReward_avg_ch1', 0,...
        'cuedReward_n_ch1', 0,...
        'cuedReward_std_ch1', 0,...
        'cuedReward_sem_ch1', 0,...
        'cuedReward_avg_ch1', 0,...
        'cuedReward_n_ch1', 0,...
        'cuedReward_std_ch1', 0,...
        'cuedReward_sem_ch1', 0 ...
    );
    


        peaks = TE.phPeakMean_cs(1).data(cuedReward);
        cComplete_summary.cuedReward_avg_ch1 = mean(peaks); % avg
        cComplete_summary.cuedReward_n_ch1 = sum(cuedReward); % n
        cComplete_summary.cuedReward_std_ch1 = std(peaks); % avg
        cComplete_summary.cuedReward_sem_ch1 = std(peaks) / sqrt(sum(cuedReward)); % SEM
% change for punish
        peaks = TE.phPeakMean_cs(1).data(cuedReward);
        cComplete_summary.cuedReward_avg_ch1 = mean(peaks); % avg
        cComplete_summary.cuedReward_n_ch1 = sum(cuedReward); % n
        cComplete_summary.cuedReward_std_ch1 = std(peaks); % avg
        cComplete_summary.cuedReward_sem_ch1 = std(peaks) / sqrt(sum(cuedReward)); % SEM        
        
        
if saveOn
    save(fullfile(savepath, ['summary_' subjectName '.mat']), 'cComplete_summary');
    disp(['*** saving: ' fullfile(savepath, ['summary_' subjectName '.mat']) ' ***']);
end
    
