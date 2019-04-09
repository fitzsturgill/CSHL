%{ 
cuedOutcome_plotAnimals
script to explore example and summary plots for the manusript
Goals:
1) Try to plot cue responses conditioned upon licking
2) Quanitfy degree of phasic vs sustained
3) lick-alignd photometry rasters
%}




%% raster plots 

DB = dbLoadExperiment('cuedOutcome');

photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
channels = [1];
climfactor = 3;


fh = [];
%%
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    if ~success
        disp('wtf');
        continue
    end    
    
    savepath = fullfile(DB.path, 'animals', animal, filesep);
    
    
%% first-lick sorted phRasters
%     ylim = [-2 8];
    xwindow = [-4 7];

    saveName = sprintf('%s_lickAligned_rasters_%s_%s', animal, photometryField, fdField);  
    fh(end + 1) =ensureFigure(saveName, 1); 
    mcLandscapeFigSetup(fh(end));
    
    lickOnsets = TE.lickLatency_cs(highValueTrials);
    lickOnsets = sort(lickOnsets);    
    subplot(1,4,1);
    eventRasterFromTE(TE, highValueTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.lickLatency_cs);
    set(gca, 'XLim', xwindow);
    title('high value'); ylabel(animal, 'Interpreter', 'none');
    
    subplot(1,4,2);
    phRasterFromTE(TE, highValueTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'sortValues', TE.lickLatency_cs, 'zeroTimes', TE.Cue, 'window', xwindow); % 'CLimFactor', CLimFactor,
    line(lickOnsets, (1:sum(highValueTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1);    
    xlabel('Time frome odor (s)');

    lickOnsets = TE.lickLatency_cs(lowValueTrials);
    lickOnsets = sort(lickOnsets);        
    subplot(1,4,3);
    eventRasterFromTE(TE, lowValueTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.lickLatency_cs);
    set(gca, 'XLim', xwindow);
    title('low value');
    
    subplot(1,4,4);
    phRasterFromTE(TE, lowValueTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'sortValues', TE.lickLatency_cs, 'zeroTimes', TE.Cue, 'window', xwindow); % 'CLimFactor', CLimFactor,             
    line(lickOnsets, (1:sum(lowValueTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1);
end
    axs = findobj(gcf, 'Type', 'axes');
    set(axs, 'FontSize', 12);
    set(axs(2:end), 'YTick', []);

h = waitbar(0, 'slowly writing pdfs');

pdf = fullfile(DB.path, 'pooled', sprintf('lickAligned_Rasters_%s_%s.pdf', photometryField, fdField));
for counter = 1:length(fh)   
    if counter == 1
        export_fig(fh(counter),pdf);  % write to pdf
    else
        export_fig(fh(counter),'-append', pdf);  % write to pdf
    end
    waitbar(counter/(length(fh)));
end

close(h);