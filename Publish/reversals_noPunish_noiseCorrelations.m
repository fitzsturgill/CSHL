% reversals_noPunish_noiseCorrelations

DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, 'pooled', filesep);
ensureDirectory(savepath);
saveOn = 1;


maxLagInSeconds = 5;
Fs = 20;
maxLag = round(maxLagInSeconds * Fs);
na = length(DB.animals);

%% initialize table to hold correlation coefficients for each animal, forget the US for now

animals = cell(length(DB.animals), 1);
[Rcue, Rcue_expFit, auROC_cue1, auROC_cue2, RCue_shift] = deal(zeros(length(DB.animals),1));
Rcue_separate = cell(length(DB.animals), 1);
Rnoise = table(animals, Rcue, RCue_shift, Rcue_expFit, auROC_cue1,...
    auROC_cue2, Rcue_separate);

Xcorrs = zeros(maxLag * 2 + 1, na);


% ensureFigure 
for counter = 1:na
    animal = DB.animals{counter};
    dbLoadAnimal(DB, animal);
    Rnoise.animals{counter} = animal;
    
    [allCue1, allUs1, allCue_expFit1, allUs_expFit1, allCue2, allUs2, allCue_expFit2, allUs_expFit2] = deal([]);
    [allCue2_shift] = deal([]);
    trialSets = [csPlusTrials & hitTrials, csMinusTrials & CRTrials];
    
    for tcounter = 1:size(trialSets, 2)
        new1 = TE.phPeakMean_optimize_cs(1).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_optimize_cs(1).data(trialSets(:, tcounter)));
        new2 = TE.phPeakMean_optimize_cs(2).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_optimize_cs(2).data(trialSets(:, tcounter)));
%         new1 = TE.phPeakMean_cs(1).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_cs(1).data(trialSets(:, tcounter)));
%         new2 = TE.phPeakMean_cs(2).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_cs(2).data(trialSets(:, tcounter)));                

        
        allCue1 = [allCue1; new1];
        allCue2 = [allCue2; new2];
        allCue2_shift = [allCue2_shift; circshift(new2, 1)];
        
        Rnoise.Rcue_separate{counter}(tcounter) = corr(new1, new2);
%         allUs1 = [allUs1; TE.phPeakMean_us(1).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_us(1).data(trialSets(:, tcounter)))];
%         allUs2 = [allUs2; TE.phPeakMean_us(2).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_us(2).data(trialSets(:, tcounter)))];

        allCue_expFit1 = [allCue_expFit1; TE.phPeakMean_cs_expFit(1).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_cs_expFit(1).data(trialSets(:, tcounter)))];
        allCue_expFit2 = [allCue_expFit2; TE.phPeakMean_cs_expFit(2).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_cs_expFit(2).data(trialSets(:, tcounter)))];
%         allUs_expFit1 = [allUs_expFit1; TE.phPeakMean_us_expFit(1).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_us_expFit(1).data(trialSets(:, tcounter)))];
%         allUs_expFit2 = [allUs_expFit2; TE.phPeakMean_us_expFit(2).data(trialSets(:, tcounter)) - nanmean(TE.phPeakMean_us_expFit(2).data(trialSets(:, tcounter)))];
    end
    Rnoise.Rcue(counter) = corr(allCue1, allCue2);
    Rnoise.Rcue_shift(counter) = corr(allCue1, allCue2_shift);
%     Rnoise.Rus(counter) = corr(allUs1, allUs2);
    Rnoise.Rcue_expFit(counter) = corr(allCue_expFit1, allCue_expFit2);
%     Rnoise.Rus_expFit(counter) = corr(allUs_expFit1, allUs_expFit2);
    Rnoise.auROC_cue1(counter) = TE.phPeakMean_optimize_cs(1).settings.optimumWindow_auROCMax;
    Rnoise.auROC_cue2(counter) = TE.phPeakMean_optimize_cs(2).settings.optimumWindow_auROCMax;    
    
    % cross correlation
    [R, lags] = avgXCorr(TE.Photometry.data(1).ZS(csPlusTrials & hitTrials, :)', TE.Photometry.data(2).ZS(csPlusTrials & hitTrials, :)', maxLag);    
    Xcorrs(:,counter) = R';            
end
lagsInSeconds = lags / Fs;

%% plot some things to check data- looks like using making sure that auROC_cue_ch1 and auROC_cue_ch1 > 0.45 selects for experiments with at least some level of noise correlations

saveName = 'Rnoise_diagnostics';
ensureFigure(saveName, 1);
colormap jet;
mcPortraitFigSetup(gcf);

subplot(2,2,1);
scatter(Rnoise.Rcue, Rnoise.Rcue_expFit, 20, 1:length(DB.animals)); addUnityLine;
xlabel('R cue'); ylabel('R cue expFit');

subplot(2,2,2);
scatter(Rnoise.Rcue.^2, Rnoise.Rcue_expFit.^2, 20, 1:length(DB.animals)); addUnityLine;
xlabel('CD cue'); ylabel('CD cue expFit');

subplot(2,2,3);
scatter(min(Rnoise.auROC_cue1, Rnoise.auROC_cue2), Rnoise.Rcue.^2, 20, 1:length(DB.animals)); 
xlabel('min cue auROC ch1&2');  ylabel('CD cue');
set(gca, 'YLim', [0 0.25]);
subplot(2,2,4);
scatter(min(Rnoise.auROC_cue1, Rnoise.auROC_cue2), Rnoise.Rcue_expFit.^2, 20, 1:length(DB.animals)); 
xlabel('min cue auROC ch1&2'); ylabel('CD cue expFit');
set(gca, 'YLim', [0 0.25]);

% subplot(3,2,5);
% scatter(min(Rnoise.auROC_cue1, Rnoise.auROC_cue2), Rnoise.Rus.^2); 
% xlabel('min cue auROC ch1&2'); ylabel('CD us');
% 
% subplot(3,2,6);
% scatter(min(Rnoise.auROC_cue1, Rnoise.auROC_cue2), Rnoise.Rus_expFit.^2); 
% xlabel('min cue auROC ch1&2'); ylabel('CD us expFit');

fig = gcf;
axs = fig.Children;
set(axs,  'XGrid', 'on', 'XMinorGrid', 'on');%, 'XLim', [0 1], 'YLim', [0 1]);

%% histogram of noise correlations for cue period
binEdges = [-1:0.1:1];
ensureFigure('noiseHist', 1); 
histogram(Rnoise.Rcue, binEdges);

cum_all = cum(Rnoise.Rcue);
cum_all_shift = cum(Rnoise.Rcue_shift);
cum_t1 = cum(cellfun(@(x) x(1), Rnoise.Rcue_separate));
cum_t2 = cum(cellfun(@(x) x(2), Rnoise.Rcue_separate));

% cumulative histogram
saveName = 'reversals_noise_cumHist';
ensureFigure(saveName, 1);
axes; plot(cum_all.sorted, cum_all.index, 'k'); hold on;
plot(cum_all_shift.sorted, cum_all_shift.index, 'Color', [0.7 0.7 0.7]);
xlabel('Rnoise');
formatFigurePublish('size', [0.8 1]);
if saveOn    
    print(gcf, '-dpdf', fullfile(savepath, filesep, 'figure', filesep, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, filesep, 'figure', filesep, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, filesep, 'figure', filesep, [saveName '.jpg']));
end
    
% plot cs+ vs cs- correlations versus each other in a scatterplot
ensureFigure('test2', 1);
scatter(cellfun(@(x) x(1), Rnoise.Rcue_separate), cellfun(@(x) x(2), Rnoise.Rcue_separate));
addUnityLine;
addOrginLines;
xlabel('Cs+'); ylabel('CS-');

%% plot cross correlations
nr = ceil(sqrt(na));
nc = ceil(na/nr);

saveName = 'reversals_xcorr_tiled';
ensureFigure(saveName, 1);
mcPortraitFigSetup(gcf);
for counter = 1:na
    subplot(nr, nc, counter);
    plot(lagsInSeconds, Xcorrs(:,counter), 'k');
    textBox(DB.animals{counter}, [], [], 8);
    set(gca, 'YLim', [-1 1]);
    xlabel('lag (s)'); ylabel('xcorr');
    addOrginLines;
end
        
        
    


%% subfunctions
function [R, lags] = avgXCorr(x, y, maxlag)
    allR = NaN(maxlag*2+1, size(x,2));
    lags = -maxlag:maxlag;

    for column = 1:size(x, 2)
        common = isfinite(x(:,column)) & isfinite(y(:,column));
        thisX = x(common,column);
        thisY = y(common,column);
        if length(thisX) ~= length(thisY)
            disp('wtf');
        end
        [theseR, lags] = xcorr(thisX, thisY, maxlag, 'coeff');
        if isempty(theseR)
            continue
        end
        allR(:,column) = theseR;
    end
    R = nanmean(allR, 2);
end
