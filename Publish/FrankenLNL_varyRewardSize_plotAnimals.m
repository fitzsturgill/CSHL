% FrankenLNL_varyRewardSize_plotAnimals


DB = dbLoadExperiment('FrankenLNL_varyRewardSize');

PhotometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
trialNumbering = 'consecutive';
CLimFactor = 3;
fhr = [];
fha = [];

for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    

    success = dbLoadAnimal(DB, animal);
    display(animal);    
    savepath = fullfile(DB.path, 'animals', animal, filesep);
    
    saveName = sprintf('Averages_%s', animal);  
    fha(end + 1) =ensureFigure(saveName, 1);
    
    trials = {trialsByType{1}, trialsByType{2}, trialsByType{3}, neutralTrials};
    window = [-0.5 3];
    
    subplot(1,3,1);
    [ha, hl] = phPlotAverageFromTE(TE, trials, 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'b', 'c', 'm', 'k'}, 'PhotometryField', PhotometryField); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    addStimulusPatch(gca, [-0.1 0.1]);
    set(gca, 'XLim', window); ylabel('Fluor ZS');
    title('Ch 1');
    t = textBox(animal, gca, [0.1 0.95]); set(t, 'Color', [1 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
        
    subplot(1,3,2);
    [ha, hl] = phPlotAverageFromTE(TE, trials, 2, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us, 'window', window, 'linespec', {'b', 'c', 'm', 'k'}, 'PhotometryField', PhotometryField); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    addStimulusPatch(gca, [-0.1 0.1]);
    set(gca, 'XLim', window); ylabel('Fluor ZS');        
    title('Ch 2');
    
    subplot(1,3,3);
    varargin = {'window', window, 'zeroTimes', TE.Us, 'window', window, 'linespec', {'b', 'c', 'm', 'k'}};
    axh = [];
    [ha, hl] = plotEventAverageFromTE(TE, trials, 'Port1In', varargin{:});  
    addStimulusPatch(gca, [-0.1 0.1]);
    ylabel('licks (1/s)'); xlabel('time from reward (s)');  set(gca, 'XLim', window);
    legend(hl, {'10', '5', '2', 'none'}, 'Box', 'off', 'Location', 'best'); 
    title('Licks');    

    % rasters
    saveName = sprintf('rasters_%s', animal); 
    fhr(end+1) = ensureFigure(saveName, 1);
    mcLandscapeFigSetup(gcf);  
    

    window = [-2 4];

    subplot(3,3,1);
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
        'zeroTimes', TE.Us, 'window', window); set(gca, 'XLim', window);         title('Licks');   
    ylabel('10uL trials');
    addStimulusPatch(gca, [-0.1 0.1]);

    subplot(3,3,2);
    phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Us);
    set(gca, 'XLim', window); title('Ch1');

    subplot(3,3,3);
    try
        phRasterFromTE(TE, trialsByType{1}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Us);
        set(gca, 'XLim', window); title('Ch2');
    end

    subplot(3,3,4);
    eventRasterFromTE(TE, trialsByType{2}, 'Port1In', 'trialNumbering', trialNumbering,...
        'zeroTimes', TE.Us, 'window', window); set(gca, 'XLim', window);         
    ylabel('5uL trials');
    addStimulusPatch(gca, [-0.1 0.1]);

    subplot(3,3,5);
    phRasterFromTE(TE, trialsByType{2}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Us);
    set(gca, 'XLim', window);

    subplot(3,3,6);
    try
        phRasterFromTE(TE, trialsByType{2}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Us);
        set(gca, 'XLim', window); 
    end

    subplot(3,3,7);
    eventRasterFromTE(TE, trialsByType{3}, 'Port1In', 'trialNumbering', trialNumbering,...
        'zeroTimes', TE.Us, 'window', window); set(gca, 'XLim', window);         
    ylabel('2uL trials'); xlabel('Time from reward (s)');
    addStimulusPatch(gca, [-0.1 0.1]);

    subplot(3,3,8);
    phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Us);
    set(gca, 'XLim', window); 

    subplot(3,3,9);
    try
        phRasterFromTE(TE, trialsByType{3}, 2, 'PhotometryField', PhotometryField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Us);
        set(gca, 'XLim', window); 
    end      
    subplot(3,3,1); t = textBox(animal, gca, [0.1 0.95]); set(t, 'Color', [1 0 0], 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left');
end

h = waitbar(0, 'slowly writing avg pdfs');
pdfavg = fullfile(DB.path, 'pooled', 'FrankenLNL_varyRewardSize_avgs_allAnimals.pdf');
for counter = 1:length(fha)    
    if counter == 1
        export_fig(fha(counter),pdfavg);  % write to pdf
    else
        export_fig(fha(counter),'-append',pdfavg);  % write to pdf
    end
    waitbar(counter/(length(fha)));
end
close(h);

h = waitbar(0, 'slowly writing raster pdfs');
pdfavg = fullfile(DB.path, 'pooled', 'FrankenLNL_varyRewardSize_rasters_allAnimals.pdf');
for counter = 1:length(fhr)    
    if counter == 1
        export_fig(fhr(counter),pdfavg);  % write to pdf
    else
        export_fig(fhr(counter),'-append',pdfavg);  % write to pdf
    end
    waitbar(counter/(length(fhr)));
end
close(h);