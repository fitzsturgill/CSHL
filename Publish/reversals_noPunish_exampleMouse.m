% reversals_noPunish_exampleMouse

% using... DC_56


DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);
smoothWindow = 1;
saveOn = 1;
animal = 'DC_56';

photometryField = 'Photometry';
fdField = 'ZS';

success = dbLoadAnimal(DB, animal); % load TE and trial lookups

%% plot rasters for an example session/reversal

    
    saveName = sprintf('example_allBehavior_csPlus_%s', animal); 
    ensureFigure(saveName, 1);
    
    csPlusExampleTrials = csPlusTrials & rewardTrials & ismember(TE.sessionIndex, [2]);
    reversals = find(diff(TE.BlockNumber(csPlusExampleTrials, :))) + 1;
    sessionChanges = find(diff(TE.sessionIndex(csPlusExampleTrials, :))) + 1;

    
    climfactor = 2;
    subplot(1,5,1);       
    CData = TE.pupil.pupDiameterNorm(csPlusExampleTrials, :) - nanmean(TE.pupil.pupDiameterNorm(csPlusExampleTrials, 1:bpX2pnt(0,20,-4)), 2);
    image(CData, 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(CData(:)) - std(CData(:), 'omitnan') * climfactor, nanmean(CData(:)) + std(CData(:), 'omitnan') * climfactor]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('pupil');  
    ylabel('trial number');
    
    
    subplot(1,5,2);
%     CData = TE.Whisk.whiskNorm(csPlusExampleTrials, :) - nanmean(TE.Whisk.whiskNorm(csPlusExampleTrials, 1:bpX2pnt(0,20,-4)), 2);
    imagesc(TE.Whisk.whiskNorm(csPlusExampleTrials, :), 'XData', [-4 7], [0 2])   
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('whisking');
    
    subplot(1,5,3);
    [~, lh] = eventRasterFromTE(TE, csPlusExampleTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7]);
    set(lh, 'LineWidth', 0.1, 'Color', [0 0 0]);
    title('licking');
    xlabel('Time from odor (s)');
    climfactor = 3;  
    subplot(1,5,4); phRasterFromTE(TE, csPlusExampleTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('ChAT'); 
    
    subplot(1,5,5); phRasterFromTE(TE, csPlusExampleTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('DAT');    
    axs = findobj(gcf, 'Type', 'axes');    
    set(axs(1:end -1), 'YTick', []);
    
    formatFigurePublish('size', [4 1.5]);
    if saveOn 
        export_fig(fullfile(savepath, saveName), '-eps');
    end
    
%% averages median split on 