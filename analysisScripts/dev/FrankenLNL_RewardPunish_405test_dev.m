%% 1/21/19 FS

saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_FrankenLNL_4odors(sessions);
%%
nTrials = sessions.SessionData.nTrials;
TE = makeTE_FrankenLNL_4odors(sessions);
TE.LED1_amp = zeros(nTrials, 1);
TE.LED2_amp = zeros(nTrials, 1);
for counter = 1:sessions.SessionData.nTrials
    TE.LED1_amp = sessions.SessionData.TrialSettings(counter).GUI.LED1_amp;
    TE.LED2_amp = sessions.SessionData.TrialSettings(counter).GUI.LED2_amp;
end
%% process photometry
% channel 1 demodulated via channel 1 reference
dFFMode = {'simple', 'simple'};
bl = [2 4];
% TE.Photometry_ch1 = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
%     'zeroField', 'Cue2', 'channels', 1, 'baseline', bl, 'downsample', 305, 'forceAmp', 1);
% 
% % channel 1 demodulated via channel 2 reference
% TE.Photometry_ch2 = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
%     'zeroField', 'Cue2', 'channels', 1, 'baseline', bl, 'downsample', 305, 'refChannels', 2, 'forceAmp', 1);

TE.Photometry_ch1 = processTrialAnalysis_Photometry_expFitConcat(sessions,...
    'zeroField', 'Cue2', 'channels', 1, 'baseline', bl, 'downsample', 305, 'forceAmp', 1);

% channel 1 demodulated via channel 2 reference
TE.Photometry_ch2 = processTrialAnalysis_Photometry_expFitConcat(sessions,...
    'zeroField', 'Cue2', 'channels', 1, 'baseline', bl, 'downsample', 305, 'refChannels', 2, 'forceAmp', 1);

duration = length(TE.Photometry_ch1.xData) / TE.Photometry_ch1.sampleRate;


%% add wheel

duration = length(TE.Photometry_ch1.xData) / TE.Photometry_ch1.sampleRate;
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', duration, 'Fs', 20, 'startField', 'PreCsRecording', 'zeroField', 'Cue2');

%%
% basepath = uigetdir;
basepath = 'Z:\SummaryAnalyses\Franken_LNL_405_dev';
sep = strfind(TE.filename{1}, '_');
if length(unique(TE.filename)) > 1
    subjectName = TE.filename{1}(1:sep(2)-1);
else
    subjectName = TE.filename{1}(1:end-4);
end
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);

%% truncate sessions
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'
TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'trialType', 1));


%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end


%% plot expFitConcat bleach fit
nPoints = numel(TE.Photometry_ch1.data(1).raw);
x = 0:nPoints-1;
fo = TE.Photometry_ch1.bleachFit.fitobject_session;
ch1_blFit = fo.a + fo.b * exp(fo.c * x);% + fo.d * exp(fo.e * x);
fo = TE.Photometry_ch2.bleachFit.fitobject_session;
ch2_blFit = fo.a + fo.b * exp(fo.c * x);% + fo.d * exp(fo.e * x);

saveName = 'expFitConcat_bleachFit';
ensureFigure(saveName, 1);
subplot(1,2,1); 
ch1data = TE.Photometry_ch1.data(1).raw';
ch1data = ch1data(:);
plot(ch1data, 'k'); hold on;
plot(ch1_blFit, 'r');
title('470nm');
subplot(1,2,2); 
ch2data = TE.Photometry_ch2.data(1).raw';
ch2data = ch2data(:);
plot(ch2data, 'k'); hold on;
plot(ch2_blFit, 'r');
title('405nm');

    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end
%% 
frankenLNL_conditions;

%% raster plots
trialNumbering = 'consecutive';
FluorDataField = 'dF';
CLimFactor = 3;
window = [-4 6];


    saveName = ['Reward_phRasters_405_test'];
    ensureFigure(saveName, 1);
    
    subplot(1,3,1);
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
        'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('cued reward');         
    addStimulusPatch(gca, [0 1]);
    subplot(1,3,2);
    phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', 'Photometry_ch1', 'FluorDataField', FluorDataField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('470nm'); set(gca, 'XLim', window);
    subplot(1,3,3);
    phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', 'Photometry_ch2', 'FluorDataField', FluorDataField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('405nm'); set(gca, 'XLim', window);
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end
    


%% averages

% -3 6
window = [-4 7];
PhotometryField = 'Photometry';
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);


    saveName = ['Averages_405test'];
    fig = ensureFigure(saveName, 1);


    axes; hold on;    
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1]), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'r'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');

    

        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1]), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'c', 'm'}, 'PhotometryField', 'Photometry_ch2'); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        set(gca, 'XLim', window); ylabel('Fluor ZS');
        
if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end
        
%%
FluorDataField = 'dF';
trialSets = [1 2 5 7];
setLabels = {'cuedReward', 'omission', 'uncuedReward', 'control'};
labels = {'baseline', 'cue', 'outcome'};
windows = [-2 0; 0 3; 3 4];
colors = [0 0 1; 0 1 0; 1 0 0];


saveName = 'distributions_trial_epochs';
h = ensureFigure(saveName, 1); colormap jet;
mcLandscapeFigSetup(h);
subplot(2,3,1);
mylines = [];
[ha, hl] = phPlotAverageFromTE(TE, trialsByType{1}, 1, 'FluorDataField', FluorDataField, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
mylines(end + 1) = hl;
[ha, hl] = phPlotAverageFromTE(TE, trialsByType{1}, 1, 'FluorDataField', FluorDataField, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'m'}, 'PhotometryField', 'Photometry_ch2'); %high value, reward
mylines(end + 1) = hl;
addStimulusPatch(gca, [0 1 0 0.8], 'odor'); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1 0 0.8], 'reinforcement');
title('cued Reward'); ylabel('dF');
legend(mylines, {'470', '405'}, 'Location', 'northwest'); legend('boxoff');

subplot(2,3,2); hold on;
[ha, hl] = phPlotAverageFromTE(TE, ~TE.reject, 1, 'FluorDataField', FluorDataField, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
[ha, hl] = phPlotAverageFromTE(TE, ~TE.reject, 1, 'FluorDataField', FluorDataField, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'m'}, 'PhotometryField', 'Photometry_ch2'); %high value, reward
addStimulusPatch(gca, [windows(1,:) 0 0.8], labels{1}, colors(1,:)); 
addStimulusPatch(gca, [windows(2,:) 0 0.8], labels{2}, colors(2,:));
addStimulusPatch(gca, [windows(3,:) 0 0.8], labels{3}, colors(3,:));

title('all trials');

        
        



ax=[];
for counter = 1:length(trialSets)
    ax(counter) = subplot(2,3,counter + 2); hold on;
    trials = trialsByType{trialSets(counter)};
    for counter2 = 1:length(labels)
        pointRange = [bpX2pnt(windows(counter2, 1), 20, -4) bpX2pnt(windows(counter2, 2), 20, -4)];
        data_470 = TE.Photometry_ch1.data(1).(FluorDataField)(trials,pointRange(1):pointRange(2));
        data_470 = data_470(:);
        data_405 = TE.Photometry_ch2.data(1).(FluorDataField)(trials,pointRange(1):pointRange(2));        
        data_405 = data_405(:);
        scatter(data_405, data_470, 20, colors(counter2, :), '.'); 
    end
    legend(labels); legend('boxoff');
    xlabel('405'); ylabel('470');
    title(setLabels(counter))
end
sameYScale(ax);
if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end


%% fit 470 to 405 using baseline period
blRange = [-2 0];
FluorDataField = 'dF'; % use dF because exponential bleaching trend across the entire session already removed from each signal
pointRange = [bpX2pnt(blRange(1), 20, -4) bpX2pnt(blRange(2), 20, -4)];
xData = TE.Photometry_ch2.data(1).(FluorDataField)(:,pointRange(1):pointRange(2));
xData = xData(:);
yData = TE.Photometry_ch1.data(1).(FluorDataField)(:,pointRange(1):pointRange(2));
yData = yData(:);

saveName = 'linearFit';
ensureFigure(saveName, 1);
scatter(xData, yData, 14, [0 0 1], '.'); hold on;
fo = fitoptions('poly1', 'Robust', 'on');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
fob = fit(xData, yData, 'poly1', fo); 
fph=plot(fob); legend off; %,'predfunc'); legend off;
set(fph, 'LineWidth', 0.5, 'Color', [1 0 0]);
textBox(sprintf('470estimated = %.2g  *  405measured + %.2g', fob.p1, fob.p2));
xlabel('405nm'); ylabel('470nm');
title('linear regression using baseline period');

% generate deltaF/F
data405 = TE.Photometry_ch2.data(1).dF;
data470 = TE.Photometry_ch1.data(1).dF;
estimated470 = data405 .* fob.p1 + fob.p2;
TE.Photometry_ch1.data(1).dF_corrected = (data470 - estimated470);% ./ estimated470; 

if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% redo rasters and averages using corrected dF
trialNumbering = 'consecutive';
CLimFactor = 3;
window = [-4 6];


    saveName = ['Reward_phRasters_dF_corrected'];
    ensureFigure(saveName, 1);
    
    subplot(1,3,1);
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
        'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('cued reward');         
    addStimulusPatch(gca, [0 1]);
    subplot(1,3,2);
    phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', 'Photometry_ch1', 'FluorDataField', 'dF_corrected', 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('corrected'); set(gca, 'XLim', window);
    subplot(1,3,3);
    phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', 'Photometry_ch1', 'FluorDataField', 'dF', 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('raw'); set(gca, 'XLim', window);
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end

    % -3 6
window = [-4 7];
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);



    saveName = ['Averages_dF_corrected'];
    fig = ensureFigure(saveName, 1);


    
    ax = subplot(1,2,1);    
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5 7]), 1, 'FluorDataField', 'dF_corrected', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'r', 'c', 'k'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');    
    set(gca, 'XLim', window); ylabel('Fluor (deltaF)'); xlabel('time from cue (s)');
    title('corrected'); 
    
    ax(2) = subplot(1,2,2);
        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5 7]), 1, 'FluorDataField', 'dF', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'r', 'c', 'k'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    
    set(gca, 'XLim', window); ylabel('Fluor (deltaF)'); xlabel('time from cue (s)');
    title('raw');
    sameYScale(ax); 
    subplot(1,2,1); addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1], '', [0.8 0.8 0.8]);addOrginLines(gca, [0 0 0]);
    subplot(1,2,2); addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1], '', [0.8 0.8 0.8]);addOrginLines(gca, [0 0 0]); 
if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end
%% compare 470 and 405
saveName = 'examples_470_vs_405';
showThese = find(Odor2Valve1Trials & rewardTrials);
ensureFigure(saveName, 1); 
whichOne = 3;
for counter = 1:4
    whichOne = counter + 10;
    subplot(2,2,counter); hold on;
    hl = zeros(2,1);
    hl(1) = plot(TE.Photometry_ch1.xData, TE.Photometry_ch1.data(1).ZS(showThese(whichOne), :)', 'b', 'LineWidth', 1); set(gca, 'XLim', [-2 4]); 
    hl(2) = plot(TE.Photometry_ch1.xData, TE.Photometry_ch2.data(1).ZS(showThese(whichOne), :)', 'm', 'LineWidth', 1); set(gca, 'XLim', [-2 4]);
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [1.9 2.1]);
    legend(hl, {'470', '405'}, 'Box', 'off', 'Location', 'northwest');
    xlabel('time from cue (s)'); ylabel('Fluor. (ZS)');
end


if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end    


% compare corrected and uncorrected
saveName = 'examples_corrected_raw';
ensureFigure(saveName, 1); 
for counter = 1:4
    whichOne = counter + 10;
    subplot(2,2,counter); hold on;
    hl = zeros(2,1);
    hl(1) = plot(TE.Photometry_ch1.xData, TE.Photometry_ch1.data(1).dF_corrected(showThese(whichOne), :)', 'b', 'LineWidth', 1); set(gca, 'XLim', [-2 4]); 
    hl(2) = plot(TE.Photometry_ch1.xData, TE.Photometry_ch1.data(1).dF(showThese(whichOne), :)', 'm', 'LineWidth', 1); set(gca, 'XLim', [-2 4]);
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [1.9 2.1]);
    legend(hl, {'corrected', 'raw'}, 'Box', 'off', 'Location', 'northwest');
    xlabel('time from cue (s)'); ylabel('Fluor. (dF)');
end


if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end    


%%
%% rasters and averages, combined
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhotometryField = 'Photometry_ch1';
FluorField = 'dF_corrected';

saveName = 'cued_reinforcement_allBehavior_plusAvgs_airPuff_405corrected';
fig = ensureFigure(saveName, 1);
mcLandscapeFigSetup(fig);

% Define the positions of different axes matrices on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    params = struct();
    matpos_title_1 = [0 0 3/8 .1];
    matpos_title_2 = [3/8 0 (1 - 3/8) .1];
    matpos_rasters = [0 .1 1 0.5];    
    matpos_avgs = [0 0.6 1 0.4];
    params.cellmargin = [.02 .02 0.05 0.05];
    params.figmargin = [0.05 0.05 0.025 0.025];
    
% add titles
params.matpos = matpos_title_1;
[ax txt] = textAxes(fig, 'Cued Reward -->', 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
set(txt, 'HorizontalAlignment', 'left', 'Position', [0 .5]);
params.matpos = matpos_title_2;
[ax txt] = textAxes(fig, sprintf('       Cued Air Puff -->    Mouse = %s', subjectName), 12);
setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);
set(txt, 'HorizontalAlignment', 'left', 'Position', [0 .5]);    

% cued Reward, cued Shock, all behavior, consecutive
trialNumbering = 'consecutive';
CLimFactor = 3;
window = [-2 5];
totalTrials = length(TE.filename);
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);


params.matpos = matpos_rasters;
nRows = 1; nColumns = 9; nAxes = 9;
hax = axesmatrix(nRows, nColumns, 1:nAxes, params, fig);

axes(hax(1)); hold on;
eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
    'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('Licks');         
addStimulusPatch(gca, [0 totalDelay]);



axes(hax(2)); hold on;
phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', PhotometryField, 'FluorDataField', FluorField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('corrected'); set(gca, 'XLim', window, 'YTick', []);



% uncorrected
axes(hax(3)); hold on;
phRasterFromTE(TE, trialsByType{1}, 1,  'PhotometryField', PhotometryField, 'FluorDataField', 'dF', 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('uncorrected'); set(gca, 'XLim', window, 'YTick', []);


try
    axes(hax(4)); hold on;
    alignedDataRaster(TE.pupil.(pupField), trialsByType{1}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
    title('Pupil diameter'); set(gca, 'XLim', window, 'YTick', []);
catch
end

axes(hax(5)); hold on;
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'FluorDataField', FluorField, 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('corrected'); set(gca, 'XLim', window);

% uncorrected
axes(hax(6)); hold on;
phRasterFromTE(TE, trialsByType{3}, 1, 'PhotometryField', PhotometryField, 'FluorDataField', 'dF', 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
title('uncorrected'); set(gca, 'XLim', window, 'YTick', []);


try
    axes(hax(7)); hold on;
    alignedDataRaster(TE.pupil.pupDiameterNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
    title('Pupil diameter'); set(gca, 'XLim', window, 'YTick', []);
catch
end

try
    axes(hax(8)); hold on;
    alignedDataRaster(TE.pupil.frameAvgNorm, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.pupil.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.pupil.frameRate(1));
    title('Eye closure'); set(gca, 'XLim', window, 'YTick', []);
catch
end

axes(hax(9)); hold on;
alignedDataRaster(TE.Wheel.data.V, trialsByType{3}, 'zeroTimes', TE.Cue2, 'window', window, 'startTimes', TE.Wheel.startTime, 'trialNumbering', trialNumbering, 'Fs', TE.Wheel.Fs, 'CLim', [0 20]);
title('Velocity'); set(gca, 'XLim', window, 'YTick', []);


for counter = 1:length(hax)
    line(hax(counter), [0 0], get(hax(counter), 'YLim'), 'Color', 'r'); line(hax(counter), [totalDelay - 0.1 totalDelay - 0.1], get(hax(counter), 'YLim'), 'Color', 'r');
end



% averages
% -3 6
params.cellmargin = [.05 .05 0.05 0.05];
params.matpos = matpos_avgs;
nRows = 1; nColumns = 4; nAxes = 4;
hax = axesmatrix(nRows, nColumns, 1:nAxes, params, fig);
% for channel = channels

% ch1, appetitive
axes(hax(1)); hold on;
[ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5]), 1,  'PhotometryField', PhotometryField, 'FluorDataField', FluorField, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k', 'c'}); %high value, reward
set(gca, 'XLim', window);  ylabel('Ch1 (dFCorrected)'); xlabel('time from cue (s)'); 
addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
lh = legend(hl, 'cued reward', 'omit', 'uncued', 'Location', 'best'); set(lh, 'Box', 'off');

% ch1, uncorrected, appetitive
axes(hax(2)); hold on;
try
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5]), 1, 'PhotometryField', PhotometryField, 'FluorDataField', 'dF', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k', 'c'}); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    set(gca, 'XLim', window);  ylabel('Ch1 (dF)');
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
catch
end

% ch1, aversive
try
    axes(hax(3)); hold on;
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([3 4 6]), 1, 'PhotometryField', PhotometryField, 'FluorDataField', FluorField, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k', 'm'}); %high value, reward
    set(gca, 'XLim', window);  ylabel('Ch1 (dFCorrected)'); xlabel('time from cue (s)');
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
    lh = legend(hl, 'cued puff', 'omit', 'uncued', 'Location', 'best'); set(lh, 'Box', 'off');
    % ch1 uncorrected, aversive
    axes(hax(4)); hold on;
catch
end
try
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([3 4 6]), 1, 'PhotometryField', PhotometryField, 'FluorDataField', 'dF', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k', 'm'}); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    set(gca, 'XLim', window);  ylabel('Ch1 (dF)');
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
catch
end






if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end


%% just plain averages


% Define the positions of different axes matrices on the figure
    % params.matpos defines position of axesmatrix [LEFT TOP WIDTH HEIGHT].    
    params = struct();
    matpos_title = [0 0 1 .1];
    matpos_avgs = [0 .1 1 0.85];    
    params.cellmargin = [.05 .05 0.05 0.05];    
    params.figmargin = [0.05 0.05 0.025 0.025];


    
    
    % -3 6
window = [-1 5];
PhotometryField = 'Photometry_ch1';
FluorField = 'dF_corrected';
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);

% rewarding subset
showTheseR = find(Odor2Valve1Trials & rewardTrials);
[~, rixR] = sort(rand(size(showTheseR)));
% aversive subset
showTheseA = find(Odor2Valve2Trials & punishTrials);
[~, rixA] = sort(rand(size(showTheseA)));


for channel = 1
    saveName = ['Averages_ch' num2str(channel)];
    fig = ensureFigure(saveName, 1);
    mcPortraitFigSetup(fig);

    params.matpos = matpos_title;
    [ax, ~] = textAxes(fig, sprintf('%s , Channel %d', subjectName, channel), 14);
    setaxesOnaxesmatrix(ax, 1, 1, 1, params, fig);    
    
    try
        ax = subplot(3,2,1);
        plot(TE.(PhotometryField).xData, TE.(PhotometryField).data(channel).(FluorField)(showTheseR(rixR(1:10)), :)', 'LineWidth', 0.25); set(gca, 'XLim', window);
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('Fluor. (ZS)'); title('rewarding');       
    catch
    end
    
    try
        ax(2) = subplot(3,2,2);
        plot(TE.(PhotometryField).xData, TE.(PhotometryField).data(channel).(FluorField)(showTheseA(rixA(1:10)), :)', 'LineWidth', 0.25); set(gca, 'XLim', window);
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('Fluor. (ZS)'); title('aversive');
    catch
    end    

    ax(3) = subplot(3,2,3);
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5]), channel, 'PhotometryField', PhotometryField, 'FluorDataField', FluorField, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k', 'c'}, 'PhotometryField', PhotometryField); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
    set(gca, 'XLim', window); ylabel('Fluor ZS');
    
    try
        ax(4) = subplot(3,2,4);
        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([3 4 6]), channel, 'PhotometryField', PhotometryField, 'FluorDataField', FluorField, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k', 'm'}, 'PhotometryField', PhotometryField); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        set(gca, 'XLim', window); ylabel('Fluor ZS');
    catch
    end
    
    try
        ax(5) = subplot(3,2,5);
        varargin = {'window', window, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'k','c'}};
        axh = [];
        [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 2 5]), 'Port1In', varargin{:});  
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('licks (1/s)'); xlabel('time from cue(s)');  set(gca, 'XLim', window);
        legend(hl, {'cued reward', 'cued ommission', 'uncued reward'}, 'Box', 'off', 'Location', 'best'); 
    catch
    end
    
    try
        ax(6) = subplot(3,2,6);
        varargin = {'window', window, 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'r', 'k','m'}};
        axh = [];
        [ha, hl] = plotEventAverageFromTE(TE, trialsByType([3 4 6]), 'Port1In', varargin{:});  
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        ylabel('licks (1/s)'); xlabel('time from cue(s)');  set(gca, 'XLim', window); legend(hl, {'cued punish', 'cued ommission', 'uncued punish'}, 'Box', 'off', 'Location', 'best');      
    catch
    end
    

    
    params.matpos = matpos_avgs;
    setaxesOnaxesmatrix(ax, 3, 2, 1:6, params, fig);     
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end    
end


