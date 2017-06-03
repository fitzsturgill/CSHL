    
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
    dFFMode{end+1} = 'simple';
    BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
    BL{end + 1} = [2 4];    
end


    

TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL);

%%
    saveName = [subjectName '_phAvgs'];  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);

    pm = [1 2];
    
    % - 6 0 4
    if ismember(1, channels)
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {filterTE(TE, 'trialType', 1), filterTE(TE, 'trialType', 2)}, 1,...
            'FluorDataField', 'dFF', 'window', [-3, 7], 'linespec', {'b', 'r'}); %high value, reward
        legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);
    end
    if ismember(2, channels)
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {filterTE(TE, 'trialType', 1), filterTE(TE, 'trialType', 2)}, 2,...
            'FluorDataField', 'dFF', 'window', [-3, 7], 'linespec', {'b', 'r'}); %high value, reward
        legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);
    end
    
    
    %%
    
    CLimFactor = 2;

reversals = find(TE.BlockChange);
for channel = channels

    saveName = [subjectName '_comboRasters_ch' num2str(channel)];
    h=ensureFigure(saveName, 1);
    mcPortraitFigSetup(h);
    

    
    subplot(1,2,1); phRasterFromTE(TE, filterTE(TE, 'trialType', 1), channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
        set(gca, 'FontSize', 14)
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines        
    
    subplot(1,2,2); phRasterFromTE(TE, filterTE(TE, 'trialType', 2), channel, 'CLimFactor', CLimFactor, 'trialNumbering', 'global');
        set(gca, 'FontSize', 14)
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines   

end
%     subplot(1,4,3); 
%     eventRasterFromTE(TE, csMinusTrials, 'Port1In', 'trialNumbering', 'global',...
%         'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'g', 'LineWidth', 2); % reversal lines    
%     title('CS-');
%     set(gca, 'XLim', [-4 7]); 
%     set(gca, 'YLim', [0 max(trialCount)]);
%     set(gca, 'FontSize', 14)