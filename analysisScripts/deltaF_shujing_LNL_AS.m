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


%% photometry averages, zscored
%     ylim = [-2 8];
    saveName = [subjectName '_phAvgs_shujing'];  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);

    pm = [2 2];
    
    % - 6 0 4
    fluorField = 'raw';
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {cuedReward}, 1,...
            'FluorDataField', fluorField, 'window', [-4, 7], 'linespec', {'g'}); %high value, reward
%         legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);
        
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {cuedReward}, 2,...
            'FluorDataField', fluorField, 'window', [-4, 7], 'linespec', {'r'}); %high value, reward
%         legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);
    
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {cuedPunish}, 1,...
            'FluorDataField', fluorField, 'window', [-4, 7], 'linespec', {'g'}); %high value, reward
%         legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);
        
        subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {cuedPunish}, 2,...
            'FluorDataField', fluorField, 'window', [-4, 7], 'linespec', {'r'}); %high value, reward
%         legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);        
        
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end    
    
    
    
    %%
%     
    TE.Photometry.data(1).dF = bsxfun(@minus, TE.Photometry.data(1).raw, TE.Photometry.data(1).blF);
    TE.Photometry.data(2).dF = bsxfun(@minus, TE.Photometry.data(2).raw, TE.Photometry.data(2).blF);
    
    %%
    
    
%     ylim = [-2 8];
    saveName = [subjectName '_phAvgs_shujing'];  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);

    pm = [1 2];
    
    % - 6 0 4
    fluorField = 'ZS';
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, trialsByType{1}, 1,...
            'FluorDataField', fluorField, 'linespec', {'g'}); %high value, reward
%         legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);
        hold on;
                
        [ha, hl] = phPlotAverageFromTE(TE, trialsByType{1}, 2,...
            'FluorDataField', fluorField, 'linespec', {'r'}); %high value, reward
%         legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);
    
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {punishTrials}, 1,...
            'FluorDataField', fluorField, 'linespec', {'g'}); %high value, reward
%         legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);
        hold on;

        [ha, hl] = phPlotAverageFromTE(TE, {punishTrials}, 2,...
            'FluorDataField', fluorField, 'linespec', {'r'}); %high value, reward
%         legend(hl, {'rew', 'pun', 'neu'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel('BF dF/F Zscored'); textBox(subjectName);%set(gca, 'YLim', ylim);        
        
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end   