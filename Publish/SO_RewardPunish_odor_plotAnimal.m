% reversals_noPunish_poolReversals

DB = dbLoadExperiment('SO_RewardPunish_odor');

photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
climfactor = 2;
window = [-6 4];

fhr = [];
fha = [];
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    
    % kludge
    if strcmp('DAT_5', animal)
        channel = 2;
    else
        channel = 1;
    end
    success = dbLoadAnimal(DB, animal);
    display(animal);
    TE.csLicks = countEventFromTE(TE, 'Port1In', [-2 0], TE.PostUsRecording);
    TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.PostUsRecording);
    TE.usZeros = cellfun(@(x,y,z) max(x(1), max(y(1), z(1))), TE.Reward, TE.Punish, TE.Omit); %'Reward', 'Punish', 'WNoise', 'Neutral'
    if ~success
        disp('wtf');
        continue
    end    
    
    savepath = fullfile(DB.path, 'animals', animal, filesep);
    
%% photometry averages, zscored
%     ylim = [-2 8];

    saveName = sprintf('Averages_%s', animal);  
    fha(end + 1) =ensureFigure(saveName, 1);
    
    subplot(2, 2, 1);
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([7 1 2 5]), channel,...
        'PhotometryField', 'Photometry', 'FluorDataField', fdField, 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'k', 'b','g','c'}); 
    legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    title('Reward'); ylabel(sprintf('BF %s', fdField)); textBox(animal, gca, [0.25 0.95], 8);%set(gca, 'YLim', ylim);
    
    subplot(2, 2, 2);
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([7 3 4 6]), channel,...
        'PhotometryField', 'Photometry', 'FluorDataField', fdField, 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'k', 'r','g','m'}); 
    legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    title('Punish'); ylabel(sprintf('BF %s', fdField)); textBox(animal, gca, [0.25 0.95], 8);%set(gca, 'YLim', ylim);
    
    subplot(2,2,3);
    varargin = {'window', [-4 0], 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'k', 'b','g','c'}};
    axh = [];
    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([7 1 2 5]), 'Port1In', varargin{:});
    legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    ylabel('licks (1/s)'); xlabel('time from reinforcement (s)'); 
    
    subplot(2,2,4);
    varargin = {'window', [-4 0], 'zeroTimes', TE.usZeros, 'window', window, 'linespec', {'k', 'r','g','m'}};
    axh = [];
    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([7 3 4 6]), 'Port1In', varargin{:});
    legend(hl, {'null', 'cued', 'omission', 'uncued'}, 'Location', 'best', 'FontSize', 8); legend('boxoff');
    ylabel('licks (1/s)'); xlabel('time from reinforcement (s)');     
    
    
    %% raster plots
    saveName = sprintf('rasters_%s', animal); 
    fhr(end+1) = ensureFigure(saveName, 1);
    mcLandscapeFigSetup(gcf);            
    trialTypes = [1 2 5 3 4 6];
    titles  = {'cued rew.', 'omit rew.', 'uncued rew.', 'cued pun.', 'omit pun.', 'uncued pun.'};
    for typeCounter = 1:length(trialTypes)
        trialType = trialTypes(typeCounter);
        sessionChanges = find(diff(TE.sessionIndex(trialsByType{trialType}, :))) + 1;
        subplot(2,6,typeCounter);
        phRasterFromTE(TE, trialsByType{trialType}, channel, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor,...
            'zeroTimes', TE.usZeros', 'window', window, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,       
        title(titles{typeCounter}); 
        
        
        subplot(2,6,typeCounter + 6);
        eventRasterFromTE(TE, trialsByType{trialType}, 'Port1In', 'trialNumbering', 'consecutive',...
            'zeroTimes', TE.usZeros, 'window', window);
        line(repmat([window(1); window(2)], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
        set(gca, 'XLim', window); 
    end
    subplot(2,6,8); xlabel('Time from Reinforcement (s)');
    subplot(2,6,1); t = textBox(animal, gca, [0.1 0.95]); set(t, 'Color', [1 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end   

h = waitbar(0, 'slowly writing avg pdfs');
pdfavg = fullfile(DB.path, 'pooled', 'SO_RewardPunish_avgs_allAnimals.pdf');
for counter = 1:length(fha)    
    if counter == 1
        export_fig(fha(counter),pdfavg);  % write to pdf
    else
        export_fig(fha(counter),'-append',pdfavg);  % write to pdf
    end
    waitbar(counter/(length(fha)*2));
end
close(h);

h = waitbar(0, 'slowly writing raster pdfs');
pdfavg = fullfile(DB.path, 'pooled', 'SO_RewardPunish_rasters_allAnimals.pdf');
for counter = 1:length(fhr)    
    if counter == 1
        export_fig(fhr(counter),pdfavg);  % write to pdf
    else
        export_fig(fhr(counter),'-append',pdfavg);  % write to pdf
    end
    waitbar(counter/(length(fhr)*2));
end
close(h);

