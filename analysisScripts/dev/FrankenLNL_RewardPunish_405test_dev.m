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
bl = [1 4];
TE.Photometry_ch1 = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Cue2', 'channels', 1, 'baseline', bl, 'downsample', 305, 'forceAmp', 1);

% channel 1 demodulated via channel 2 reference
TE.Photometry_ch2 = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
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

%% 
frankenLNL_conditions;

%% raster plots
PhotometryField = 'Photometry';
trialNumbering = 'consecutive';
CLimFactor = 3;
window = [-4 6];


    saveName = ['Reward_phRasters_405_test'];
    ensureFigure(saveName, 1);
    
    subplot(1,3,1);
    eventRasterFromTE(TE, trialsByType{1}, 'Port1In', 'trialNumbering', trialNumbering,...
        'zeroTimes', TE.Cue2, 'window', window); set(gca, 'XLim', window);         xlabel('time from cue (s)');  title('cued reward');         
    addStimulusPatch(gca, [0 1]);
    subplot(1,3,2);
    phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', 'Photometry_ch1', 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('470nm'); set(gca, 'XLim', window);
    subplot(1,3,3);
    phRasterFromTE(TE, trialsByType{1}, 1, 'PhotometryField', 'Photometry_ch2', 'CLimFactor', CLimFactor, 'trialNumbering', trialNumbering, 'window', window, 'zeroTimes', TE.Cue2);
    title('405nm'); set(gca, 'XLim', window);
    if saveOn
        saveas(gcf, fullfile(savepath, saveName), 'fig');
        saveas(gcf, fullfile(savepath, saveName), 'jpeg');
    end
    


%% averages

% -3 6
window = [-1 5];
PhotometryField = 'Photometry';
totalDelay = TE.Trace2{1}(2) - TE.Cue2{1}(1);


    saveName = ['Averages_405test'];
    fig = ensureFigure(saveName, 1);
    mcPortraitFigSetup(fig);

axes; hold on;    
    [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1]), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'g'}, 'PhotometryField', 'Photometry_ch1'); %high value, reward
%     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');

    

        [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1]), 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue2, 'window', window, 'linespec', {'m'}, 'PhotometryField', 'Photometry_ch2'); %high value, reward
    %     legend(hl, {'pRhigh', 'pRmedium', 'pRlow'}, 'Location', 'northwest', 'FontSize', 12); legend('boxoff');
        addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [totalDelay - 0.1 totalDelay + 0.1]);
        set(gca, 'XLim', window); ylabel('Fluor ZS');
        
%%

data_470 = TE.Photometry_ch1.data(1).ZS(:, 20:end - 20);
data_470 = data_470(:);
data_405 = TE.Photometry_ch2.data(1).ZS(:, 20:end - 20);
data_405 = data_405(:);

ensureFigure('test_scatter', 1);
scatter(data_405, data_470, '.'); xlabel('405'); ylabel('470');

%% do robust regression, look at residuals

