
DB = dbLoadExperiment('FrankenLNL_RewardPunish');

photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
trialNumbering = 'consecutive';
climfactor = 3;
fhr = [];
minRewardLickRate = 2; % at least n Hz licking (during us window for reward trials)
minCueLickRate = 1; % at least n Hz licking (during cs window for reward cued reward trials)
xwindow = [-2 5];
savePath = fullfile(DB.path, 'pooled');
%%


xlim = [-6 4];
for counter = 1:length(DB.animals)
% for counter = 2
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
    display(animal);
    
    trials = rewardTrials & Odor2Valve1Trials & (TE.licks_us.rate > minRewardLickRate) & (TE.licks_cs.rate > minCueLickRate) & ismember(TE.BlockNumber, [2 3])...
        & (round(cellfun(@(x) diff(x), TE.Trace2)) == max((round(cellfun(@(x) diff(x), TE.Trace2))))); % exclude shock or extinction days
    lickOnsets = TE.lickLatency_cs(trials);
    lickOnsets = sort(lickOnsets);    
    lickZeros = TE.lickLatency_cs + cellfun(@(x) x(1), TE.Cue2);
    
    saveName = sprintf('alignedReward_rasters_%s', animal);  
    fhr(end + 1) = ensureFigure(saveName, 1); 
    mcPortraitFigSetup(gcf);
    for ch = 1:2
        subplot(2,2,1 + (ch-1) * 2); hold on;
        phRasterFromTE(TE, trials, ch, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'zeroTimes', TE.Cue2, 'window', xwindow, 'sortValues', TE.lickLatency_cs); % 'CLimFactor', CLimFactor,
        line(lickOnsets, (1:sum(trials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--');    
%         set(gca, 'YTickLabel', {});
        xlabel('Time frome odor (s)');
        ylabel(['Ch ' num2str(ch)]);
        if ch == 1
            t = textBox(animal); set(t, 'Color', [1 1 1], 'FontSize', 16, 'FontWeight', 'bold');
        end



        subplot(2,2,2 + (ch-1) * 2); hold on;
        phRasterFromTE(TE, trials, ch, 'zeroTimes', lickZeros, 'window', xwindow,...
            'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'showSessionBreaks', 0, 'sortValues', TE.lickLatency_cs); % 'CLimFactor', CLimFactor,
        line(-1 * lickOnsets, (1:sum(trials))', 'Parent', gca, 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');    
        % scatter(-1 * TE.lickLatency_cs(highValueTrials & rewardTrials), 1:sum(highValueTrials & rewardTrials), markerSize, [1 1 1], 'filled');
        xlabel('Time frome lick (s)');
        ylabel(['Ch ' num2str(ch)]);
%         set(gca, 'YTickLabel', {});        
    end
end
%
h = waitbar(0, 'slowly writing pdfs');
pdfPlus = fullfile(DB.path, 'pooled', 'FrankenLNL_alignedReward_rasters.pdf');
for counter = 1:length(fhr)    
    if counter == 1
        export_fig(fhr(counter),pdfPlus);  % write to pdf
    else
        export_fig(fhr(counter),'-append',pdfPlus);  % write to pdf
    end
    waitbar(counter/length(fhr));
end
close(h);