% load('TE.mat');
%plot(nanmean(TE.Photometry.data.ZS,1));
%path = 'C:\Users\tcare\Documents\GitHub\CSHLPartners\Deconvolution\';
% path = 'C:\Users\tcare\OneDrive\Documents\GitHub\CSHLPartners\Deconvolution\';
path = 'C:\Users\Adam\Documents\Repos\CSHLPartners\Deconvolution';
cd(path);
load('TE.mat');
cuedOutcome_Conditions;

kernel = phAverageFromTE(TE, punishTrials & uncuedTrials, 1, 'window', [0 2], 'FluorDataField', 'ZS');
data =  TE.Photometry.data.ZS;
h = kernel.Avg; 
out = deconv_Fourier(data, h, 0.2);

TE.Photometry.data.Deconv = out;

% 
% full = figure('units','normalized','outerposition',[0 0 1 1]);
% hold on;
% avged = nanmean(data);
% plot(avged);% / max(avged));
% plot(out);
% legend("Before", "After");
% axis([0 180 -2 2]);
% title("Final Output");

%% averages raw

    % plot photometry averages
    h=ensureFigure('Photometry_Averages', 1); 
    mcLandscapeFigSetup(h);

    pm = [3 2];
    
    % - 6 0 4
    subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 3 7]), 1, 'FluorDataField', 'ZS'); %high value, reward
    legend(hl, {'hival, rew', 'hival, omit', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('high value'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); textBox(TE.filename{1}(1:7));

    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 6 8]), 1, 'FluorDataField', 'ZS'); % low value, punish
    legend(hl, {'loval, pun', 'loval, omit', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('low value'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1, 'FluorDataField', 'ZS'); % reward, varying degrees of expectation
    legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('reward all'); ylabel('Z Score'); xlabel('time from reinforcement (s)');     

    subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 2 8]), 1, 'FluorDataField', 'ZS'); % punishment, varying degrees of expectation
    legend(hl, {'loval, pun', 'hival, pun', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('punish all'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 5, 'FontSize', 12, 'LineWidth', 1); [ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
        'window', [-6 0], 'linespec', {'m', 'g'}, 'FluorDataField', 'ZS'); hold on;
    
    subplot(pm(1), pm(2), 5); [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, omitTrials}, 1,...
        'window', [0 4], 'linespec', {'b', 'r', 'k'}, 'FluorDataField', 'ZS');
    hl = [hla hl];
    legend(hl, {'loval', 'hival', 'rew', 'pun', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Balazs'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 
    
    subplot(pm(1), pm(2), 6, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([3 6 9]), 1, 'FluorDataField', 'ZS'); % reward, varying degrees of expectation
    legend(hl, {'hival, neutral', 'loval, neutral', 'neutral'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('neutral all'); ylabel('Z Score'); xlabel('time from reinforcement (s)');    
    
if saveOn    
    saveas(gcf, fullfile(savepath, 'phAverages.fig'));
    saveas(gcf, fullfile(savepath, 'phAverages.jpg'));
end

%%
    % plot photometry averages
    h=ensureFigure('Photometry_Averages_deconv', 1); 
    mcLandscapeFigSetup(h);

    pm = [3 2];
    
    % - 6 0 4
    subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 3 7]), 1, 'FluorDataField', 'Deconv'); %high value, reward
    legend(hl, {'hival, rew', 'hival, omit', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('high value'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); textBox(TE.filename{1}(1:7));

    subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 6 8]), 1, 'FluorDataField', 'Deconv'); % low value, punish
    legend(hl, {'loval, pun', 'loval, omit', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('low value'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1, 'FluorDataField', 'Deconv'); % reward, varying degrees of expectation
    legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('reward all'); ylabel('Z Score'); xlabel('time from reinforcement (s)');     

    subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([5 2 8]), 1, 'FluorDataField', 'Deconv'); % punishment, varying degrees of expectation
    legend(hl, {'loval, pun', 'hival, pun', 'pun'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('punish all'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 

    subplot(pm(1), pm(2), 5, 'FontSize', 12, 'LineWidth', 1); [ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
        'window', [-6 0], 'linespec', {'m', 'g'}, 'FluorDataField', 'Deconv'); hold on;
    
    subplot(pm(1), pm(2), 5); [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials, omitTrials}, 1,...
        'window', [0 4], 'linespec', {'b', 'r', 'k'}, 'FluorDataField', 'Deconv');
    hl = [hla hl];
    legend(hl, {'loval', 'hival', 'rew', 'pun', 'omit'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('Balazs'); ylabel('Z Score'); xlabel('time from reinforcement (s)'); 
    
    subplot(pm(1), pm(2), 6, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([3 6 9]), 1, 'FluorDataField', 'Deconv'); % reward, varying degrees of expectation
    legend(hl, {'hival, neutral', 'loval, neutral', 'neutral'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
    title('neutral all'); ylabel('Z Score'); xlabel('time from reinforcement (s)');    
    
if saveOn    
    saveas(gcf, fullfile(savepath, 'phAverages.fig'));
    saveas(gcf, fullfile(savepath, 'phAverages.jpg'));
end