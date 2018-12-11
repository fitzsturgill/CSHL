% initial_learning_plotAnimal


DB = dbLoadExperiment('initial_learning');

photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
channels = [1 2];
climfactor = 2;

fhr = [];
fha = [];
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    if ~success
        disp('wtf');
        continue
    end    
    
    savepath = fullfile(DB.path, 'animals', animal, filesep);
    
    %% raster plots
    CLimFactor = 2;
   
    saveName = sprintf('initial_learning_rasters_%s', animal);  
    fhr(end + 1) =ensureFigure(saveName, 1);

    subplot(1,4,1); 
    eventRasterFromTE(TE, csPlusTrials & (day1Trials | day2PlusTrials) & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroTimes', TE.Cue, 'window', [-4 7]);%'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    title('cued'); ylabel('trial number');
    set(gca, 'YLim', [0 sum(csPlusTrials & (day1Trials | day2PlusTrials) & rewardTrials)]);
    set(gca, 'XLim', [-4 7]); 
    set(gca)


    subplot(1,4,2);
    phRasterFromTE(TE, csPlusTrials & (day1Trials | day2PlusTrials) & rewardTrials, 1, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'zeroTimes', TE.Cue, 'window', [-4 7]);
    set(gca);
    xlabel('Time from cue (s)');

    subplot(1,4,3);
    phRasterFromTE(TE, csPlusTrials & (day1Trials | day2PlusTrials) & rewardTrials, 2, 'CLimFactor', CLimFactor, 'trialNumbering', 'consecutive', 'zeroTimes', TE.Cue, 'window', [-4 7]);
    set(gca);    
    xlabel('Time from cue (s)');
    subplot(1,4,4); 
    textAxes(gca, {'cued reward', animal}, 10);
end



h = waitbar(0, 'slowly writing pdfs');

pdfRasters = fullfile(DB.path, 'pooled', sprintf('initial_learning_rasters.pdf'));
for counter = 1:length(fhr)    
    if counter == 1
        export_fig(fhr(counter),pdfRasters);  % write to pdf
    else
        export_fig(fhr(counter),'-append',pdfRasters);  % write to pdf
    end
    waitbar(counter/(length(fhr)));
end
close(h);