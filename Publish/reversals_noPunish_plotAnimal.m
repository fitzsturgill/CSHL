% reversals_noPunish_poolReversals

DB = dbLoadExperiment('reversals_noPunish_publish');


saveOn = 1;
channels = [1 2];

fhp = [];
fhm = [];
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    if ~success
        disp('wtf');
        continue
    end    
    
    savepath = fullfile(DB.path, 'animals', animal, filesep);
    
%% photometry averages, zscored
%     ylim = [-2 8];
    fdField = 'ZS';
    saveName = sprintf('%s_phAvgs_%s', animal, fdField);  
    h=ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(h);

    pm = [2 2];
    
    % - 6 0 4
    if ismember(1, channels)
        subplot(pm(1), pm(2), 1, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, neutralTrials, uncuedReward}, 1,...
            'FluorDataField', fdField, 'window', [1, 7], 'linespec', {'b','k','c'}); %high value, reward
        legend(hl, {'rew', 'neutral', 'uncued rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('Reinforcement'); ylabel(sprintf('BF %s', fdField)); textBox(animal);%set(gca, 'YLim', ylim);
    end
    
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 3, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, neutralTrials, uncuedReward}, 2,...
            'FluorDataField', fdField, 'window', [1, 7], 'linespec', {'b','k','c'}); %high value, reward
        legend(hl, {'rew', 'neutral', 'uncued rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        ylabel(sprintf('VTA %s', fdField)); xlabel('time from cue (s)'); %set(gca, 'YLim', ylim);
    end
    
    % - 6 0 4
    if ismember(1, channels)    
        subplot(pm(1), pm(2), 2, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 1,...
        'FluorDataField', fdField, 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        title('CS+, outcomes'); set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);
    end
    if ismember(2, channels)    
        subplot(pm(1), pm(2), 4, 'FontSize', 12, 'LineWidth', 1); 
        [ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 2,...
            'FluorDataField', fdField, 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
        legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
        xlabel('time from cue (s)');     set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);
    end
    
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end        
       
%%  pupil averages    
    saveName = sprintf('Pupil_avgs_%s', animal); 
    ensureFigure(saveName, 1); 
    % condition on behavior and cue condition
    subplot(2,2,1);[~, hl] = plotPupilAverageFromTE(TE, {uncuedTrials, csPlusTrials & hitTrials, csMinusTrials & CRTrials}, 'window', [-2 6]);
    legend(hl, {'uncued', 'cs+, hit', 'cs-, CR'}, 'Box', 'off', 'Location', 'northwest');
    title('Cue condition'); xlabel('time from cue (s)');
    set(gca, 'XLim', [-2 6]);
    % condition on behavior only for csPlus
    subplot(2,2,2);[~, hl] = plotPupilAverageFromTE(TE, {hitTrials & csPlusTrials, missTrials & csPlusTrials}, 'window', [-2 6]);
    set(gca, 'XLim', [-2 6]);
    title('Cs+ by behavioral response'); xlabel('time from cue (s)');
    legend(hl, {'cs+, hit', 'cs+, miss'}, 'Box', 'off', 'Location', 'northwest');
    % reward, expected vs unexpected
    subplot(2,2,3);[~, hl] = plotPupilAverageFromTE(TE, {csPlusTrials & hitTrials & rewardTrials, (csMinusTrials & CRTrials & rewardTrials) | uncuedReward}, 'window', [-2 6]);
    set(gca, 'XLim', [-2 6]);
    legend(hl, {'expected reward', 'unexpected reward'}, 'Box', 'off', 'Location', 'northwest');
    title('Reward by cue and behavior'); xlabel('time from cue (s)');
    formatFigure('scaleFactor', 4)
    if saveOn
        saveas(gcf, fullfile(savepath, [saveName '.fig']));
        saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    end
    
%% all behavior CsPlus
    fdField = 'ZS';
    saveName = sprintf('allBehavior_csPlus_%s', animal); 
    fhp(end+1) = ensureFigure(saveName, 1);
    
    reversals = find(diff(TE.BlockNumber(csPlusTrials, :))) + 1;
    sessionChanges = find(diff(TE.sessionIndex(csPlusTrials, :))) + 1;

    
    
    subplot(1,6,1);    
    image(TE.Wheel.data.V(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 3, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 3]); 
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1);
    t = textBox(animal); set(t, 'Color', [1 1 1], 'FontSize', 16, 'FontWeight', 'bold');
    title('Velocity');
    ylabel('trial number');

    subplot(1,6,2);
    
    try
        image(TE.pupil.pupDiameterNorm(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
        set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
        line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
        colormap('parula');  
        title('Pupil Diameter');    
    catch
        disp('wtf');        
    end
    title('pupil');
    
    
    subplot(1,6,3);
    imagesc(TE.Whisk.whiskNorm(csPlusTrials, :), 'XData', [-4 7], [0 2])
    
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('whisking');
    
    subplot(1,6,4);
    eventRasterFromTE(TE, csPlusTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7]);
    title('licking');
    
    
    subplot(1,6,5); phRasterFromTE(TE, csPlusTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 1, 'FluorDataField', fdField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,6,6); phRasterFromTE(TE, csPlusTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 1, 'FluorDataField', fdField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('DAT');    
    axs = findobj(gcf, 'Type', 'axes');

    set(axs, 'FontSize', 16);
    set(axs(2:end), 'YTick', []);
    set(gcf, 'Position', [1 1 1920 1004]);
    saveas(gcf, fullfile(savepath, saveName), 'fig'); 
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    
%% all behavior CsMinus
    saveName = sprintf('allBehavior_csMinus_%s', animal);     
    fhm(end + 1) = ensureFigure(saveName, 1);
    
    reversals = find(diff(TE.BlockNumber(csMinusTrials, :))) + 1;
    sessionChanges = find(diff(TE.sessionIndex(csMinusTrials, :))) + 1;
    
    subplot(1,6,1);    
    image(TE.Wheel.data.V(csMinusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines      
    title('Velocity');
    ylabel('trial number');
    t = textBox(animal); set(t, 'Color', [1 1 1], 'FontSize', 16, 'FontWeight', 'bold');

    subplot(1,6,2);    
    try
        image(TE.pupil.pupDiameterNorm(csMinusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
        set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
        line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
        colormap('parula');  
        title('Pupil Diameter');    
    catch
    end
    title('pupil');
    
    
    subplot(1,6,3);
    imagesc(TE.Whisk.whiskNorm(csMinusTrials, :), 'XData', [-4 7], [0 2])
    
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('whisking');
    
    subplot(1,6,4);
    eventRasterFromTE(TE, csMinusTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7]);
    title('licking');
    
    
    subplot(1,6,5); phRasterFromTE(TE, csMinusTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    
    subplot(1,6,6); phRasterFromTE(TE, csMinusTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('DAT');    
    axs = findobj(gcf, 'Type', 'axes');

    set(axs, 'FontSize', 16);
    set(axs(2:end), 'YTick', []);
    set(gcf, 'Position', [1 1 1920 1004]);
    saveas(gcf, fullfile(savepath, 'allBehavior_whisk_csMinus'), 'fig'); 
    saveas(gcf, fullfile(savepath, 'allBehavior_whisk_csMinus'), 'jpeg');
end

% save pdf versions of allBehavior

h = waitbar(0, 'slowly writing pdfs');

pdfPlus = fullfile(DB.path, 'pooled', 'allBehavior_csPlus.pdf');
for counter = 1:length(fhp)    
    if counter == 1
        export_fig(fhp(counter),pdfPlus);  % write to pdf
    else
        export_fig(fhp(counter),'-append',pdfPlus);  % write to pdf
    end
    waitbar(counter/length(fhp));
end

pdfMinus = fullfile(DB.path, 'pooled', 'allBehavior_csMinus.pdf');
for counter = 1:length(fhm)    
    if counter == 1
        export_fig(fhm(counter),pdfMinus);  % write to pdf
    else
        export_fig(fhm(counter),'-append',pdfMinus);  % write to pdf
    end
    waitbar((counter + length(fhm))/length(fhm) * 2);
end
close(h);
    