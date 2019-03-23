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
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 5]), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'r'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');

    

        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 5]), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'c', 'm'}, 'PhotometryField', 'Photometry_ch2'); %high value, reward
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
labels = {'post', 'cue', 'outcome'};
windows = [4 6.9; 0 3; 3 4];
colors = [1 0 0; 0 1 0; 0 0 1];

saveName = 'test_scatter';
ensureFigure(saveName, 1); colormap jet;
ax=[];
for counter = 1:length(trialSets)
    ax(counter) = subplot(2,2,counter); hold on;
    trials = trialsByType{trialSets(counter)};
    for counter2 = 1:length(labels)
        pointRange = [bpX2pnt(windows(counter2, 1), 20, -4) bpX2pnt(windows(counter2, 2), 20, -4)];
        data_470 = TE.Photometry_ch1.data(1).(FluorDataField)(trials,pointRange(1):pointRange(2));
        data_470 = data_470(:);
        data_405 = TE.Photometry_ch2.data(1).(FluorDataField)(trials,pointRange(1):pointRange(2));        
        data_405 = data_405(:);
        scatter(data_405, data_470, 15, colors(counter2, :), '.'); 
    end
    legend(labels);
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
scatter(xData, yData, 14, [0 0 0], '.'); hold on;
fo = fitoptions('poly1', 'Robust', 'on');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
fob = fit(xData, yData, 'poly1', fo); 
fph=plot(fob); legend off; %,'predfunc'); legend off;
set(fph, 'LineWidth', 0.5, 'Color', [1 0 0]);

% generate deltaF/F
data405 = TE.Photometry_ch2.data(1).dF;
data470 = TE.Photometry_ch1.data(1).dF;
estimated470 = data405 .* fob.p1 + fob.p2;
TE.Photometry_ch1.data(1).dF_corrected = (data470 - estimated470);% ./ estimated470; 

%% redo rasters and averages using corrected dFF
trialNumbering = 'consecutive';
FluorDataField = 'dFF_corrected';
CLimFactor = 3;
window = [-4 6];


    saveName = ['Reward_phRasters_dFF_corrected'];
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
PhotometryField = 'Photometry';
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);



    saveName = ['Averages_dF_corrected'];
    fig = ensureFigure(saveName, 1);


    
    ax = subplot(1,2,1);    
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5 7]), 1, 'FluorDataField', 'dF_corrected', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'r', 'c', 'k'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
    set(gca, 'XLim', window); ylabel('Fluor (deltaF)'); xlabel('time from cue (s)');
    title('corrected');
    
    ax(2) = subplot(1,2,2);
        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 2 5 7]), 1, 'FluorDataField', 'dF', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'b', 'r', 'c', 'k'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
    addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
    set(gca, 'XLim', window); ylabel('Fluor (deltaF)'); xlabel('time from cue (s)');
    title('raw');
    sameYScale(ax);
if saveOn
    saveas(gcf, fullfile(savepath, saveName), 'fig');
    saveas(gcf, fullfile(savepath, saveName), 'jpeg');
end
%% synchrony between left and right BLA
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


