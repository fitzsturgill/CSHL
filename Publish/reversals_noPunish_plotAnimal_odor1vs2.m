% reversals_noPunish_poolReversals

DB = dbLoadExperiment('reversals_noPunish_publish');

photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
channels = [1 2];
climfactor = 2;

fh = [];
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    if ~success
        disp('wtf');
        continue
    end    
    
    savepath = fullfile(DB.path, 'animals', animal, filesep);
    

    
%% all behavior odor1 vs 2
    saveName = sprintf('allBehavior_Odor1vs2_%s_%s_%s', animal, photometryField, fdField); 
    fh(end+1) = ensureFigure(saveName, 1);
    
    reversals1 = find(diff(TE.BlockNumber(Odor1Trials, :))) + 1;
    sessionChanges1 = find(diff(TE.sessionIndex(Odor1Trials, :))) + 1;
    reversals2 = find(diff(TE.BlockNumber(Odor2Trials, :))) + 1;
    sessionChanges2 = find(diff(TE.sessionIndex(Odor2Trials, :))) + 1;
    
    

    
    subplot(1,6,1);
    eventRasterFromTE(TE, Odor1Trials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges1)), [sessionChanges1'; sessionChanges1'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals1)), [reversals1'; reversals1'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7]);
    title('Odor 1');
    t = textBox(animal); set(t, 'Color', [0 0 0], 'FontSize', 16, 'FontWeight', 'bold');
    
    
    subplot(1,6,2); phRasterFromTE(TE, Odor1Trials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals1)), [reversals1'; reversals1'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    set(gca, 'XLim', [-4 7], 'YTick', []);
    title('ChAT'); xlabel('Time frome odor (s)');
    
    subplot(1,6,3); phRasterFromTE(TE, Odor1Trials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals1)), [reversals1'; reversals1'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    set(gca, 'XLim', [-4 7], 'YTick', []);
    title('DAT');    
    
    subplot(1,6,4);
    eventRasterFromTE(TE, Odor2Trials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges2)), [sessionChanges2'; sessionChanges2'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals2)), [reversals2'; reversals2'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7]);
    title('Odor 2');
    
    
    subplot(1,6,5); phRasterFromTE(TE, Odor2Trials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals2)), [reversals2'; reversals2'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    set(gca, 'XLim', [-4 7], 'YTick', []);
    title('ChAT'); xlabel('Time frome odor (s)');
    
    subplot(1,6,6); phRasterFromTE(TE, Odor2Trials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals2)), [reversals2'; reversals2'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    set(gca, 'XLim', [-4 7], 'YTick', []);
    title('DAT');    
    
    axs = findobj(gcf, 'Type', 'axes');

    set(axs, 'FontSize', 16);    
    set(gcf, 'Position', [1 1 1920 1004]);
    saveas(gcf, fullfile(savepath, saveName), 'fig'); 
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');    
end

% save pdf versions of all odors

h = waitbar(0, 'slowly writing all odors pdfs');
pdfBoth = fullfile(DB.path, 'pooled', sprintf('allBehavior_Odor1vs2_%s_%s.pdf', photometryField, fdField));
for counter = 1:length(fh)    
    if counter == 1
        export_fig(fh(counter),pdfBoth);  % write to pdf
    else
        export_fig(fh(counter),'-append',pdfBoth);  % write to pdf
    end
    waitbar(counter/length(fh));
end
close(h);
