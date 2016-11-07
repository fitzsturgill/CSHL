% cuedOutcome_Odor_Complete_analysisScript_addPupil
%%

% TE = addPupilometryToTE_special(TE, 'duration', [11; 11; 10; 10; 10; 10], 'normMode', 'byTrial');
TE = addPupilometryToTE_special(TE, 'normMode', 'byTrial');
%%
phasicWindow = [-2.7, -2.2];
sustainedWindow = [-1.5, 0];
TE.phasicAvg= bpCalcPeak_dFF(TE.Photometry, 1, phasicWindow, TE.Us, 'method', 'mean');
TE.sustainedAvg= bpCalcPeak_dFF(TE.Photometry, 1, sustainedWindow, TE.Us, 'method', 'mean');

TE.phasicLicks = countEventFromTE(TE, 'Port1In', phasicWindow, TE.Us);
TE.sustainedLicks = countEventFromTE(TE, 'Port1In', sustainedWindow, TE.Us);

pupLag = 0.3; % pupil dilation lags cholinergic signal (Nelson and Mooney)
TE.phasicPup = bpCalcPeak_Pupil(TE.pupil, phasicWindow + pupLag, TE.Us);
TE.sustainedPup = bpCalcPeak_Pupil(TE.pupil, sustainedWindow + pupLag, TE.Us);
TE.pupBaseline = bpCalcPeak_Pupil(TE.pupil, [1 4]);

TE.cueCondition = ismember(TE.trialType, 1:3) * 2 + ismember(TE.trialType, 4:6) * 1;

%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
% basepath = uigetdir;
% subjectName = TE.filename{1}(1:7);
% disp(subjectName);
% savepath = fullfile(basepath, subjectName, '_pupil');
% ensureDirectory(savepath);

%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end
%%
% ensureFigure('pupTest', 1);
% 
% 
% 
% x1 = TE.sustainedPup.data(lowValueTrials);
% y1 = TE.phasicAvg.data(lowValueTrials);
% x = x1(~isnan(x1) & ~isnan(y1));
% y = y1(~isnan(x1) & ~isnan(y1));
% f = fit(x, y, 'poly1');
% plot(f, x, y, 'predfunc');
% xlabel('pupil diameter');
% ylabel('dFF sustained');

%% 
%%
% ensureFigure('pupTest', 1);
% % scatter(TE.sustainedPup.data(highValueTrials), TE.phasicAvg.data(highValueTrials), 'b'); hold on;
% % scatter(TE.sustainedPup.data(lowValueTrials), TE.phasicAvg.data(lowValueTrials), 'r'); hold on;
% % scatter(TE.sustainedPup.data(omitTrials), TE.phasicAvg.data(omitTrials), 'g'); hold on;
% x1 = TE.sustainedPup.data(lowValueTrials);
% y1 = TE.sustainedAvg.data(lowValueTrials);
% x = x1(~isnan(x1) & ~isnan(y1));
% y = y1(~isnan(x1) & ~isnan(y1));
% f = fit(x, y, 'poly1');
% plot(f, x, y, 'predfunc');
% xlabel('pupil diameter');
% ylabel('dFF sustained');

%% lets try and regress out the cue and then fit the residuals to the pupil measurements


trialsToTest = validTrials;
sustainedPup = TE.sustainedPup.data(trialsToTest);
sustainedLicks = TE.sustainedLicks.rate(trialsToTest);
sustainedPh = TE.sustainedAvg.data(trialsToTest);
cueCondition = categorical(TE.cueCondition(trialsToTest))'; % assert categorical data type
baselinePup = TE.pupBaseline.data(trialsToTest);

regressData = table(cueCondition, sustainedPh); 
% default is for last table entry to be the response variable
mdl = fitlm(regressData, 'linear');
%
residuals = mdl.Residuals.Raw;

h = ensureFigure('sustainedCue_modelFig', 1);
mcLandscapeFigSetup(h);
mdl_pup = fitlm(table(sustainedPup, residuals), 'linear');
% mdl_pup = fitlm(table(baselinePup, residuals), 'linear');
mdl_licks = fitlm(table(sustainedLicks, residuals), 'linear');
mdl_licks2 = fitlm(table(sustainedLicks, residuals), 'linear', 'Exclude', mdl_licks.Diagnostics.CooksDistance > 1);
subplot(2,2,1);
plot(mdl);
subplot(2,2,2); 
plot(mdl_pup);
subplot(2,2,3);
plot(mdl_licks2);

if saveOn
    saveas(gcf, fullfile(savepath, 'sustainedCue_modelFig.fig'));
    saveas(gcf, fullfile(savepath, 'sustainedCue_modelFig.jpg'));
end

%% do the same for phasic component

phasicPup = TE.phasicPup.data(validTrials);
phasicLicks = TE.phasicLicks.rate(validTrials);
phasicPh = TE.phasicAvg.data(validTrials);
cueCondition = categorical(TE.cueCondition(validTrials))'; % assert categorical data type

regressData = table(cueCondition, phasicPh); 
% default is for last table entry to be the response variable
mdlPhasic = fitlm(regressData, 'linear');

residuals = mdlPhasic.Residuals.Raw;

ensureFigure('Residual_fits_phasic', 1);
mdlPhasic_pup = fitlm(table(phasicPup, residuals), 'linear');
mdlPhasic_licks = fitlm(table(phasicLicks, residuals), 'linear');
subplot(2,2,1);
plot(mdlPhasic);
subplot(2,2,2); 
plot(mdlPhasic_pup);
subplot(2,2,3);
plot(mdlPhasic_licks);

if saveOn
    saveas(gcf, fullfile(savepath, 'phasicCue_modelFig.fig'));
    saveas(gcf, fullfile(savepath, 'phasicCue_modelFig.jpg'));
end




%%
%%
pupilZero = bpX2pnt(7, 60, 0);
% phZero = bpX2pnt(6, 20, 0);
ensureFigure('pupilAverages_cue', 1);
subplot(1,2,1);
phPlotAverageFromTE(TE, {highValueTrials, lowValueTrials, uncuedTrials}, 1, 'linespec', {'b', 'r', 'k'}, 'window', [-7 0]); hold on;
set(gca, 'XLim', [-7 0], 'XGrid', 'on');set(gca, 'FontSize', 12, 'TickDir', 'out');
title('ChAT dF/F (hi vs lo value cue)');
subplot(1,2,2); 
plot(TE.pupil.xData(1:pupilZero), nanmean(TE.pupil.pupDiameterNorm(highValueTrials, (1:pupilZero))), 'b'); % 
hold on;
plot(TE.pupil.xData(1:pupilZero), nanmean(TE.pupil.pupDiameterNorm(lowValueTrials, (1:pupilZero))), 'r'); % 
plot(TE.pupil.xData(1:pupilZero), nanmean(TE.pupil.pupDiameterNorm(uncuedTrials, (1:pupilZero))), 'k'); % 
set(gca, 'XLim', [-7 0], 'XGrid', 'on');set(gca, 'FontSize', 12, 'TickDir', 'out');
xlabel('time from reinforcement (s)');
ylabel('Pupil Diameter norm (by session)');
title('Pupil (hi vs lo value cue)');

if saveOn
    saveas(gcf, fullfile(savepath, 'pupilAverages_cue.fig'));
    saveas(gcf, fullfile(savepath, 'pupilAverages_cue.jpg'));
end

%% wtf number of frames lower in certain sessions
nSessions = max(TE.sessionIndex);
% ensureFigure('test', 1);
% for counter = 1:nSessions
%     subplot(3,2, counter);
%     histogram(cellfun(@(x) x.nFrames, TE.pupil.settings(TE.sessionIndex == counter)));
%     title(num2str(counter));    
% end

ensureFigure('test2', 1);
xdata = linspace(-7, 4, 660);
for counter = 1:nSessions
    subplot(3,2,counter);
    plot(xdata, nanmean(TE.pupil.eyeAreaNorm(lowValueTrials & TE.sessionIndex == counter, :))); 
    set(gca, 'XLim', [-7 4]);
    title(num2str(counter));
end

