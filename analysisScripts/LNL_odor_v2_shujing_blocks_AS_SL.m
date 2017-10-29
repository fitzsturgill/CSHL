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
    dFFMode{end+1} = 'expFit';
    BL{end + 1} = [1 4];    
end
  
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL);
%% extrat peaks either with dFF or ZS
%%
phField = 'dFF';
%%
phField = 'ZS';

%% extract peak trial dFF responses to cues and reinforcement and lick counts
% zero is defined as time of cue- see call to
% processTrialAnalysis_Photometry2
nTrials = length(TE.filename);

csWindow = [0 1];
baselineWindow = [-1 0];
TE.csLicks = countEventFromTE(TE, 'Port1In', csWindow, TE.Cue);

usWindow = [0 1];


% line below extracts time when the us occurs (accounting for the various
% outcome state possibilities)
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'
TE.Us = usZeros;
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros); %wider window for counting US licks than photometry US response

for channel = channels
    TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, BL{channel}, [], 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMedium_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');
    TE.phPeakMedium_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');
    TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.9, 'phField', 'ZS');
    TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, usWindow, usZeros, 'method', 'percentile', 'percentile', 0.9, 'phField', 'ZS');
    TE.phPeakPercentile_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, baselineWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.9, 'phField', 'ZS');
%     TE.phPeakPercentile_cs_baselined(channel) = TE.phPeakPercentile_cs(channel).data - TE.phPeakPercentile_baseline(channel).data
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
    
%     rewardTrials = filterTE(TE, 'trialType', [1], 'reject', 0); 
%     punishTrials = filterTE(TE, 'trialType', [3], 'reject', 0); 
% this part is wrong for blocks 2

    rewardTrials_SL = filterTE(TE, 'ReinforcementOutcome', 'Reward', 'reject', 0);
    punishTrials_SL = filterTE(TE, 'ReinforcementOutcome', 'Punish', 'reject', 0);
    
    cuedReward = filterTE(TE, 'OdorValveIndex', 1, 'ReinforcementOutcome', 'Reward', 'reject', 0);
    uncuedReward = filterTE(TE, 'OdorValveIndex', 0, 'ReinforcementOutcome', 'Reward', 'reject', 0);
    rewardOmission = filterTE(TE, 'OdorValveIndex', 1, 'ReinforcementOutcome', 'Neutral', 'reject', 0);
    cuedPunish = filterTE(TE, 'OdorValveIndex', 2, 'ReinforcementOutcome', 'Punish', 'reject', 0);
    punishOmission = filterTE(TE, 'OdorValveIndex', 2, 'ReinforcementOutcome', 'Neutral', 'reject', 0);
    uncuedPunish = filterTE(TE, 'OdorValveIndex', 0, 'ReinforcementOutcome', 'Punish', 'reject', 0);
    
    firstNTrials = 20;
    
    nSessions = max(TE.sessionIndex);
    Odor1FirstNTrials = [];
    Odor2FirstNTrials = [];
    for counter = 1:nSessions
        trialsThisSession = find(Odor1Trials & (TE.sessionIndex == counter));
        Odor1FirstNTrials = [Odor1FirstNTrials; trialsThisSession(1:min(firstNTrials, length(trialsThisSession)))];
        trialsThisSession = find(Odor2Trials & (TE.sessionIndex == counter));
        Odor2FirstNTrials = [Odor2FirstNTrials; trialsThisSession(1:min(firstNTrials, length(trialsThisSession)))];        
    end
        
%     Odor2TrialsIx = find(Odor2Trials);
%     
%     % odor 2
%     theseSessionIndices = TE.sessionIndex(Odor2Trials);
%     sessionOffsets = find(diff(theseSessionIndices)) + 1; % contains the indices of the first trial in a new session WITHIN Odor1Trials
%     Odor2FirstNTrials = Odor2TrialsIx(1:20);
%     for counter = 1:length(sessionOffsets)
%         offset = sessionOffsets(counter);
%         Odor2FirstNTrials = [Odor2FirstNTrials; Odor2TrialsIx(offset:offset+firstNTrials - 1)];
%     end
    
    % NOTE INCONSISTENCY ACROSS BLOCK 1 AND 2 OF TRIAL TYPES-  IN BLOCK 2
    % TRIAL TYPE 3 = CUED PUNISHMENT WHEREAS IN BLOCK 1 = UNCUED REWARD
    trialTypes = 1:6; % max 6 trial types when using shujing_blocks, block #2

    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end

    trialCount = [1:length(TE.filename)]';
    
    %% lick and photometry rasters
    
    
    CLimFactor = 2; % scale each image +/-  2 * global baseline standard deviation

for channel = channels
    saveName = [subjectName '_rasters_ch' num2str(channel)];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);
    
    subplot(1,3,1); % lick raster
    eventRasterFromTE(TE, cuedReward, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('Cued Reward'); ylabel('trial number');
    set(gca, 'XLim', [-4 6]); 
%     set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)        
    
    % photometry raster
    subplot(1,3,2); phRasterFromTE(TE, cuedReward, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        set(gca, 'XLim', [-4 6]); set(gca, 'FontSize', 14)

    subplot(1,3,3); 
    if sum(cuedPunish)
        phRasterFromTE(TE, cuedPunish, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        title('Cued Punish');
        set(gca, 'XLim', [-4 6]); set(gca, 'FontSize', 14)  
    end
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
        'binWidth', 0.5, 'window', [-4, 6], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'linespec', {'b', 'k'}};
    axh = [];
    subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, {cuedReward, rewardOmission}, 'Port1In', varargin{:});
    legend(hl, {'cued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Licks'); ylabel('licks (s)'); xlabel('time from cue (s)'); 
    
    
    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
    % channel provided as the 3rd argument
    [ha, hl] = phPlotAverageFromTE(TE, {cuedReward, rewardOmission}, 1,...
        'FluorDataField', 'ZS', 'window', [-4, 6], 'linespec', {'b', 'k'}); %high value, reward
    legend(hl, {'cued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Ch1'); ylabel('Ch1 ZS');
                        
    if ismember(2, channels)
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        % channel provided as the 3rd argument
        [ha, hl] = phPlotAverageFromTE(TE, {cuedReward, rewardOmission}, 2,...
            'FluorDataField', 'ZS', 'window', [-4, 6], 'linespec', {'b', 'k'}); %high value, reward
        legend(hl, {'cued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Ch2'); ylabel('Ch2 ZS');               
    end
    
    
%     varargin = {'trialNumbering', 'consecutive',...
%         'binWidth', 0.5, 'window', [-4, 6], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'linespec', {'b', 'c'}};
%     axh = [];
%     subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, {cuedRewardFirst, uncuedReward}, 'Port1In', varargin{:});
%     legend(hl, {'cued', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
%     title('Licks'); ylabel('licks (s)'); xlabel('time from cue (s)'); 
%     
%     
%     subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
%     % channel provided as the 3rd argument
%     [ha, hl] = phPlotAverageFromTE(TE, {cuedRewardFirst, uncuedReward}, 1,...
%         'FluorDataField', 'ZS', 'window', [-4, 6], 'linespec', {'b', 'c'}); %high value, reward
%     legend(hl, {'cued', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
%     title('Ch1'); ylabel('Ch1 dF/F');
%                         
%     if ismember(2, channels)
%         subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
%         % channel provided as the 3rd argument
%         [ha, hl] = phPlotAverageFromTE(TE, {cuedRewardFirst, uncuedReward}, 2,...
%             'FluorDataField', 'ZS', 'window', [-4, 6], 'linespec', {'b', 'c'}); %high value, reward
%         legend(hl, {'cued', 'uncued'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
%         title('Ch2'); ylabel('Ch2 dF/F');               
%     end
    
    
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
        'binWidth', 0.5, 'window', [-4, 6], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'linespec', {'r', 'k'}};
    axh = [];
    subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, {cuedPunish, punishOmission}, 'Port1In', varargin{:});
    legend(hl, {'cued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Licks'); ylabel('licks (s)'); xlabel('time from cue (s)'); 
    
    
    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
    % channel provided as the 3rd argument
    [ha, hl] = phPlotAverageFromTE(TE, {cuedPunish, punishOmission}, 1,...
        'FluorDataField', 'ZS', 'window', [-4, 6], 'linespec', {'r', 'k'}); %high value, reward
    legend(hl, {'cued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Ch1'); ylabel('Ch1 ZS');
                        
    if ismember(2, channels)
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        % channel provided as the 3rd argument
        [ha, hl] = phPlotAverageFromTE(TE, {cuedPunish, punishOmission}, 2,...
            'FluorDataField', 'ZS', 'window', [-4, 6], 'linespec', {'r', 'k'}); %high value, reward
        legend(hl, {'cued', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Ch2'); ylabel('Ch2 ZS');             
    end
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));   % for illustrator, uncomment        
    end  

%% Averages, Reward & Punish
    saveName = [subjectName '_reward & punish Avgs'];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);

    pm = [2 1]; % subplot matrix
    
    varargin = {'trialNumbering', 'consecutive',...
        'binWidth', 0.5, 'window', [-4, 6], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'linespec', {'r', 'k'}};
    axh = [];
    
      subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
    % channel provided as the 3rd argument
    [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials_SL, punishTrials_SL}, 1,...
        'FluorDataField', 'ZS', 'window', [-4, 6], 'linespec', {'r', 'k'}); %high value, reward
    legend(hl, {'reward', 'punish'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Ch1'); ylabel('Ch1 ZS');
   
                        
    if ismember(2, channels)
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        % channel provided as the 3rd argument
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials_SL, punishTrials_SL}, 2,...
            'FluorDataField', 'ZS', 'window', [-4, 6],  'linespec', {'r', 'k'}); %high value, reward
        legend(hl, {'reward', 'punish'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Ch2'); ylabel('Ch2 ZS');             
    end
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));   % for illustrator, uncomment        
    end  
    
    
    %% Average, Odor first n trials
    saveName = [subjectName '_odor_First' num2str(firstNTrials) 'Trial_Avgs'];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);

    pm = [3 1]; % subplot matrix
    
    varargin = {'trialNumbering', 'consecutive',...
        'binWidth', 0.5, 'window', [-4, 6], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'linespec', {'b', 'c'}};
    axh = [];
    subplot(pm(1), pm(2), 1); [ha, hl] = plotEventAverageFromTE(TE, {Odor1FirstNTrials}, 'Port1In', varargin{:});
    legend(hl, {'odor1FirstNTrials'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Licks'); ylabel('licks (s)'); xlabel('time from cue (s)'); 
    
    
    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
    % channel provided as the 3rd argument
    [ha, hl] = phPlotAverageFromTE(TE, {Odor1FirstNTrials, Odor2FirstNTrials}, 1,...
        'FluorDataField', 'ZS', 'window', [-4, 6], 'linespec', {'b', 'c'}); %high value, reward
    legend(hl, {'Odor1FirstNTrials', 'Odor2FirstNTrials'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Ch1'); ylabel('Ch1 ZS');
    
    if ismember(2, channels)
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
            [ha, hl] = phPlotAverageFromTE(TE, {Odor1FirstNTrials, Odor2FirstNTrials}, 2,...
            'FluorDataField', 'ZS', 'window', [-4, 6], 'linespec', {'b', 'c'}); %high value, reward
        legend(hl, {'Odor1FirstNTrials', 'Odor2FirstNTrials'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Ch2'); ylabel('Ch2 ZS');               
    end
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));           
    end
    
    
%% Raster, Odor first n trials  
  for channel = channels
    CLimFactor = 2;
    saveName = [subjectName '_odor_First' num2str(firstNTrials) 'Trial_rasters' num2str(channel)];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);
    subplot(1,3,1); % lick raster
    eventRasterFromTE(TE, Odor1FirstNTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('Odor1FirstNTrials'); ylabel('trial number');
    set(gca, 'XLim', [-4 6]); 
%     set(gca, 'YLim', [0 max(trialCount)]);
    set(gca, 'FontSize', 14)        
    
    % photometry raster
    subplot(1,3,2); phRasterFromTE(TE, Odor1FirstNTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
        set(gca, 'XLim', [-4 6]); set(gca, 'FontSize', 14)

    subplot(1,3,3); phRasterFromTE(TE, Odor2FirstNTrials, channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive');
        title('Odor2FirstNTrials');
        set(gca, 'XLim', [-4 6]); set(gca, 'FontSize', 14)  
   
  
  
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
              
    end  
    
  end
  
    %% save odor 1 session peak and lick data
    odor1Data.cs_peakPercentile_ch1 = TE.phPeakPercentile_cs(1).data(Odor1FirstNTrials);
    odor1Data.csLicks  =  TE.csLicks.rate(Odor1Trials);
    
    sjbla_71_odor1 = odor1Data

    save(fullfile(savepath, [subjectName '_odor1_First' num2str(firstNTrials) 'Trials_novelty.mat']), 'sjbla_71_odor1');
    disp(['*** saving: ' fullfile(savepath, [subjectName '_odor1_First' num2str(firstNTrials) 'Trials_novelty.mat']) ' ***']);
    
    saveName = [subjectName '_odor1_First' num2str(firstNTrials) 'Trials_novelty'];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);

    %ensureFigure('novelty', 1);
    plot(odor1Data.cs_peakPercentile_ch1);

%    save(fullfile(savepath, ['odor1_' subjectName '.mat']), 'odor1Data');
%    disp(['*** saving: ' fullfile(savepath, ['odor1_' subjectName '.mat']) ' ***']);
   
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
              
    end

    %% save odor 2 session peak and lick data
    odor2Data.cs_peakPercentile_ch1 = TE.phPeakPercentile_cs(1).data(Odor2FirstNTrials);
    
    sjbla_71_odor2 = odor2Data
    
    save(fullfile(savepath, [subjectName '_odor2_First' num2str(firstNTrials) 'Trials_novelty.mat']), 'sjbla_71_odor2');
    disp(['*** saving: ' fullfile(savepath, [subjectName '_odor2_First' num2str(firstNTrials) 'Trials_novelty.mat']) ' ***']);
    
    saveName = [subjectName '_odor2_First' num2str(firstNTrials) 'Trials_novelty'];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);

    %ensureFigure('novelty', 1);
    plot(odor2Data.cs_peakPercentile_ch1);

%    save(fullfile(savepath, ['odor1_' subjectName '.mat']), 'odor1Data');
%    disp(['*** saving: ' fullfile(savepath, ['odor1_' subjectName '.mat']) ' ***']);
   
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
        saveas(gcf, fullfile(savepath, [saveName '.epsc']));   % for illustrator, uncomment        
    end
    
    
 %% save reward session peak
%     rewardpeakData.us_peakMean_ch1 = TE.phPeakMean_us(1).data(cuedReward);
%     
%     sjbla_71_rewardpeak = rewardpeakData
%     
%     save(fullfile(savepath, [subjectName '_rewardpeak'  '.mat']), 'sjbla_71_rewardpeak');
%     disp(['*** saving: ' fullfile(savepath, [subjectName '_rewardpeak'  '.mat']) ' ***']);
%    
%     %plot 
%     saveName = [subjectName '_rewardpeak'];
%     h=ensureFigure(saveName, 1);
%     mcPortraitFigSetup(h);
% 
%     %ensureFigure('rewardpeak', 1);
%     plot(rewardpeakData.us_peakMean_ch1);
%     
%     if saveOn
%         saveas(gcf, fullfile(savepath, [saveName '.fig']));
%         saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
%         saveas(gcf, fullfile(savepath, [saveName '.epsc']));   % for illustrator, uncomment        
%     end
%% summary statistics
   cComplete_summary = struct(...
       'cuedReward_n_ch1', 0,...
       'cuedReward_mean_avg_ch1', 0,...       
       'cuedReward_mean_std_ch1', 0,...
       'cuedReward_mean_sem_ch1', 0,...
       'cuedReward_medium_avg_ch1', 0,...       
       'cuedReward_medium_std_ch1', 0,...
       'cuedReward_medium_sem_ch1', 0,...
       'cuedReward_perc_avg_ch1', 0,...       
       'cuedReward_perc_std_ch1', 0,...
       'cuedReward_perc_sem_ch1', 0,...
       'cuedRewardBaselined_avg_ch1', 0,...
       'cuedRewardBaselined_std_ch1', 0,...
       'cuedRewardBaselined_sem_ch1', 0,...       
       'cuedPunish_n_ch1', 0,...
       'cuedPunish_mean_avg_ch1', 0,...       
       'cuedPunish_mean_std_ch1', 0,...
       'cuedPunish_mean_sem_ch1', 0,...
       'cuedPunish_medium_avg_ch1', 0,...       
       'cuedPunish_medium_std_ch1', 0,...
       'cuedPunish_medium_sem_ch1', 0,...
       'cuedPunish_perc_avg_ch1', 0,...       
       'cuedPunish_perc_std_ch1', 0,...
       'cuedPunish_perc_sem_ch1', 0,...
       'cuedPunishBaselined_avg_ch1', 0,...
       'cuedPunishBaselined_std_ch1', 0,...
       'cuedPunishBaselined_sem_ch1', 0,...     
       'odor1_n_ch1', 0,...
       'odor1_mean_avg_ch1', 0,...       
       'odor1_mean_std_ch1', 0,...
       'odor1_mean_sem_ch1', 0,...
       'odor1_medium_avg_ch1', 0,...       
       'odor1_medium_std_ch1', 0,...
       'odor1_medium_sem_ch1', 0,...
       'odor1_perc_avg_ch1', 0,...       
       'odor1_perc_std_ch1', 0,...
       'odor1_perc_sem_ch1', 0,...
       'odor1Baselined_avg_ch1', 0,...
       'odor1Baselined_std_ch1', 0,...
       'odor1Baselined_sem_ch1', 0,...
       'odor2_n_ch1', 0,...
       'odor2_mean_avg_ch1', 0,...       
       'odor2_mean_std_ch1', 0,...
       'odor2_mean_sem_ch1', 0,...
       'odor2_medium_avg_ch1', 0,...       
       'odor2_medium_std_ch1', 0,...
       'odor2_medium_sem_ch1', 0,...
       'odor2_perc_avg_ch1', 0,...       
       'odor2_perc_std_ch1', 0,...
       'odor2_perc_sem_ch1', 0,...
       'odor2Baselined_avg_ch1', 0,...
       'odor2Baselined_std_ch1', 0,...
       'odor2Baselined_sem_ch1', 0);
     % mean   
       mean_reward = TE.phPeakMean_us(1).data(cuedReward);
       cComplete_summary.cuedReward_n_ch1 = sum(cuedReward); % n
       cComplete_summary.cuedReward_mean_avg_ch1 = mean(mean_reward); % avg
       cComplete_summary.cuedReward_mean_std_ch1 = std(mean_reward); % std
       cComplete_summary.cuedReward_mean_sem_ch1 = std(mean_reward) / sqrt(sum(cuedReward)); % SEM
     
     % medium  
       medium_reward = TE.phPeakMedium_us(1).data(cuedReward);
       cComplete_summary.cuedReward_medium_avg_ch1 = mean(medium_reward); % avg
       cComplete_summary.cuedReward_medium_std_ch1 = std(medium_reward); % avg
       cComplete_summary.cuedReward_medium_sem_ch1 = std(medium_reward) / sqrt(sum(cuedReward)); % SEM       
             
     % percentile  0.9
       Percentiles_reward = TE.phPeakPercentile_us(1).data(cuedReward);
       cComplete_summary.cuedReward_n_ch1 = sum(cuedReward); % n
       cComplete_summary.cuedReward_perc_avg_ch1 = mean(Percentiles_reward); % avg
       cComplete_summary.cuedReward_perc_std_ch1 = std(Percentiles_reward); % avg
       cComplete_summary.cuedReward_perc_sem_ch1 = std(Percentiles_reward) / sqrt(sum(cuedReward)); % SEM
     
     % percentile cuedReward-baseline   
       baseline_reward = TE.phPeakPercentile_baseline(1).data(cuedReward);
       cComplete_summary.cuedRewardBaselined_avg_ch1 = mean(Percentiles_reward - baseline_reward); % avg
       cComplete_summary.cuedRewardBaselined_std_ch1 = std(Percentiles_reward - baseline_reward); % avg
       cComplete_summary.cuedRewardBaselined_sem_ch1 = std(Percentiles_reward - baseline_reward) / sqrt(sum(cuedReward)); % SEM
               
    % change for punish
     % mean   
       mean_punish = TE.phPeakMean_us(1).data(cuedPunish);
       cComplete_summary.cuedPunish_n_ch1 = sum(cuedPunish); % n
       cComplete_summary.cuedPunish_mean_avg_ch1 = mean(mean_punish); % avg
       cComplete_summary.cuedPunish_mean_std_ch1 = std(mean_punish); % std
       cComplete_summary.cuedPunish_mean_sem_ch1 = std(mean_punish) / sqrt(sum(cuedPunish)); % SEM
     
     % medium  
       medium_punish = TE.phPeakMedium_us(1).data(cuedPunish);
       cComplete_summary.cuedPunish_medium_avg_ch1 = mean(medium_punish); % avg
       cComplete_summary.cuedPunish_medium_std_ch1 = std(medium_punish); % avg
       cComplete_summary.cuedPunish_medium_sem_ch1 = std(medium_punish) / sqrt(sum(cuedPunish)); % SEM       
             
     % percentile  0.9
       Percentiles_punish = TE.phPeakPercentile_us(1).data(cuedPunish);
       cComplete_summary.cuedPunish_n_ch1 = sum(cuedPunish); % n
       cComplete_summary.cuedPunish_perc_avg_ch1 = mean(Percentiles_punish); % avg
       cComplete_summary.cuedPunish_perc_std_ch1 = std(Percentiles_punish); % avg
       cComplete_summary.cuedPunish_perc_sem_ch1 = std(Percentiles_punish) / sqrt(sum(cuedPunish)); % SEM
     
     % percentile cuedPunish-baseline   
       baseline_punish = TE.phPeakPercentile_baseline(1).data(cuedPunish);
       cComplete_summary.cuedPunishBaselined_avg_ch1 = mean(Percentiles_punish - baseline_punish); % avg
       cComplete_summary.cuedPunishBaselined_std_ch1 = std(Percentiles_punish - baseline_punish); % avg
       cComplete_summary.cuedPunishBaselined_sem_ch1 = std(Percentiles_punish - baseline_punish) / sqrt(sum(cuedPunish)); % SEM
 
    % change for odor1
     % mean   
       mean_odor1 = TE.phPeakMean_cs(1).data(Odor1FirstNTrials);
       cComplete_summary.odor1_n_ch1 = numel(Odor1FirstNTrials); % n
       cComplete_summary.odor1_mean_avg_ch1 = mean(mean_odor1); % avg
       cComplete_summary.odor1_mean_std_ch1 = std(mean_odor1); % std
       cComplete_summary.odor1_mean_sem_ch1 = std(mean_odor1) / sqrt(numel(Odor1FirstNTrials)); % SEM
     
     % medium  
       medium_odor1 = TE.phPeakMedium_cs(1).data(Odor1FirstNTrials);
       cComplete_summary.odor1_medium_avg_ch1 = mean(medium_odor1); % avg
       cComplete_summary.odor1_medium_std_ch1 = std(medium_odor1); % avg
       cComplete_summary.odor1_medium_sem_ch1 = std(medium_odor1) / sqrt(numel(Odor1FirstNTrials)); % SEM       
             
     % percentile  0.9
       Percentiles_odor1 = TE.phPeakPercentile_cs(1).data(Odor1FirstNTrials);
       cComplete_summary.odor1_perc_avg_ch1 = mean(Percentiles_odor1); % avg
       cComplete_summary.odor1_perc_std_ch1 = std(Percentiles_odor1); % avg
       cComplete_summary.odor1_perc_sem_ch1 = std(Percentiles_odor1) / sqrt(numel(Odor1FirstNTrials)); % SEM
     
     % percentile odor1-baseline   
       baseline_odor1 = TE.phPeakPercentile_baseline(1).data(Odor1FirstNTrials);
       cComplete_summary.odor1Baselined_avg_ch1 = mean(Percentiles_odor1 - baseline_odor1); % avg
       cComplete_summary.odor1Baselined_std_ch1 = std(Percentiles_odor1 - baseline_odor1); % avg
       cComplete_summary.odor1Baselined_sem_ch1 = std(Percentiles_odor1 - baseline_odor1) / sqrt(numel(Odor1FirstNTrials)); % SEM
       
     % change for odor2
     % mean   
       mean_odor2 = TE.phPeakMean_cs(1).data(Odor2FirstNTrials);
       cComplete_summary.odor2_n_ch1 = numel(Odor2FirstNTrials); % n
       cComplete_summary.odor2_mean_avg_ch1 = mean(mean_odor2); % avg
       cComplete_summary.odor2_mean_std_ch1 = std(mean_odor2); % std
       cComplete_summary.odor2_mean_sem_ch1 = std(mean_odor2) / sqrt(numel(Odor2FirstNTrials)); % SEM
     
     % medium  
       medium_odor2 = TE.phPeakMedium_cs(1).data(Odor2FirstNTrials);
       cComplete_summary.odor2_medium_avg_ch1 = mean(medium_odor2); % avg
       cComplete_summary.odor2_medium_std_ch1 = std(medium_odor2); % avg
       cComplete_summary.odor2_medium_sem_ch1 = std(medium_odor2) / sqrt(numel(Odor2FirstNTrials)); % SEM       
             
     % percentile  0.9
       Percentiles_odor2 = TE.phPeakPercentile_cs(1).data(Odor2FirstNTrials);
       cComplete_summary.odor2_perc_avg_ch1 = mean(Percentiles_odor2); % avg
       cComplete_summary.odor2_perc_std_ch1 = std(Percentiles_odor2); % avg
       cComplete_summary.odor2_perc_sem_ch1 = std(Percentiles_odor2) / sqrt(numel(Odor2FirstNTrials)); % SEM
     
     % percentile odor2-baseline   
       baseline_odor2 = TE.phPeakPercentile_baseline(1).data(Odor2FirstNTrials);
       cComplete_summary.odor2Baselined_avg_ch1 = mean(Percentiles_odor2 - baseline_odor2); % avg
       cComplete_summary.odor2Baselined_std_ch1 = std(Percentiles_odor2 - baseline_odor2); % avg
       cComplete_summary.odor2Baselined_sem_ch1 = std(Percentiles_odor2 - baseline_odor2) / sqrt(numel(Odor2FirstNTrials)); % SEM
       
%      % mean 
%        peaks_punish = TE.phPeakMean_us(1).data(cuedPunish);
%        cComplete_summary.cuedPunish_avg_ch1 = mean(peaks_punish); % avg
%        cComplete_summary.cuedPunish_n_ch1 = sum(cuedPunish); % n
%        cComplete_summary.cuedPunish_std_ch1 = std(peaks_punish); % avg
%        cComplete_summary.cuedPunish_sem_ch1 = std(peaks_punish) / sqrt(sum(cuedPunish)); % SEM    
%      % percentile   
%        peakPercentiles_punish = TE.phPeakPercentile_us(1).data(cuedPunish);
%        cComplete_summary2.cuedPunish_avg_ch1 = mean(peakPercentiles_punish); % avg
%        cComplete_summary2.cuedPunish_n_ch1 = sum(cuedPunish); % n
%        cComplete_summary2.cuedPunish_std_ch1 = std(peakPercentiles_punish); % avg
%        cComplete_summary2.cuedPunish_sem_ch1 = std(peakPercentiles_punish) / sqrt(sum(cuedPunish)); % SEM
%     
%        % percentile cuedPunish-baseline   
%        peakPercentiles_baseline = TE.phPeakPercentile_baseline(1).data(cuedPunish);
%        cComplete_summary2.cuedPunishBaselined_avg_ch1 = mean(peakPercentiles_punish - peakPercentiles_baseline); % avg
%        cComplete_summary2.cuedPunishBaselined_n_ch1 = sum(cuedPunish); % n
%        cComplete_summary2.cuedPunishBaselined_std_ch1 = std(peakPercentiles_punish - peakPercentiles_baseline); % avg
%        cComplete_summary2.cuedPunishBaselined_sem_ch1 = std(peakPercentiles_punish - peakPercentiles_baseline) / sqrt(sum(cuedPunish)); % SEM       
       
if saveOn
   save(fullfile(savepath, ['summary_' subjectName '_' phField '.mat']), 'cComplete_summary');
   disp(['*** saving: ' fullfile(savepath, ['summary_' subjectName '.mat']) ' ***']);
end
