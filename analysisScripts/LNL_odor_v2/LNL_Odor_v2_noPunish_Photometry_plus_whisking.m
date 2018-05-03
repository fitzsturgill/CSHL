% 4/10/17  Analysis script for pavlovian reversals using LickNoLick_Odor_V2
% protocol

saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_LNL_odor_V2(sessions);
%%
% assume that photometry channels are consistent across sessions, bleach
% fit dFF for GCaMP6f (ch1) and simple dFF for jRGECO1a (ch2)
channels=[]; dFFMode = {}; BL = {}; 
if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'expFit';
%     dFFMode{end+1} = 'simple';    
    BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
    BL{end + 1} = [0 4];    
end

%% baseline by trial
 TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL,...
    'tau', 2);   
%% baseline expfit
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', {'expFit', 'expFit'}, 'zeroField', 'Cue', 'channels', channels, 'baseline', BL,...
    'tau', 2);
% TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'expFitBegin', 0.1,...
%     'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels, 'baseline', BL, 'downsample', 305);
%%
% if you are reloading TE do this:
channels = [1 2];

%% extract peak trial dFF responses to cues and reinforcement and lick counts
% csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);
csWindow = zeros(nTrials, 2);
csWindow(:,2) = cellfun(@(x,y,z) max(x(end), y(end)) - z(1), TE.AnswerLick, TE.AnswerNoLick, TE.Cue); 
% max 1) to select AnswerNoLick time stampfor no lick trials (unused state contains NaN)
% 2) to select AnswerLick time stamp for lick trials (AnswerLick follows
% AnswerNoLick state)

ch1CsWindow = [0.25 1];
ch2CsWindow = [0.25 1];

TE.csLicks = countEventFromTE(TE, 'Port1In', [0 3], TE.Cue);

usWindow = [0 0.5];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'


for channel = channels
    TE.phPeakMean_baseline(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [1 4], [], 'method', 'mean', 'phField', 'ZS');
    TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'mean', 'phField', 'ZS');
    if channel == 1
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch1CsWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');
    elseif channel == 2
        TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'mean', 'phField', 'ZS');
        TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, ch2CsWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');        
        TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.75], usZeros, 'method', 'percentile', 'percentile', 0.5, 'phField', 'ZS');  
    end
end

% TE.wheel_cs = mean(TE.Wheel.data.V(:,bpX2pnt(ch1CsWindow(1) + pupLag, 20, -4):bpX2pnt(ch1CsWindow(2) + pupLag, 20, -4)), 2);


TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);
%%
TE = addPupilometryToTE(TE, 'duration', 11, 'zeroField', 'Cue',  'frameRate', 60, 'frameRateNew', 20);
pupLag = 0.3;

TE.pupilBaseline = mean(TE.pupil.pupDiameterNorm(:,bpX2pnt(-3, 20, -4):bpX2pnt(0, 20, -4)), 2);
% ch1CsWindow = [0.25 1];
TE.pupil_cs = mean(TE.pupil.pupDiameterNorm(:,bpX2pnt(ch1CsWindow(1) + pupLag, 20, -4):bpX2pnt(ch1CsWindow(2) + pupLag, 20, -4)), 2);
%%

TE.whisk = addWhiskingToTE(TE, 'folderPrefix', 'WhiskDiff_');

%%
phasicWindow = [-2.7, -2.2];
sustainedWindow = [-1.5, 0];
TE.phasicAvg= bpCalcPeak_dFF(TE.Photometry, 1, phasicWindow, TE.Us, 'method', 'mean');
TE.sustainedAvg = bpCalcPeak_dFF(TE.Photometry, 1, sustainedWindow, TE.Us, 'method', 'mean');
TE.phasicLicks = countEventFromTE(TE, 'Port1In', phasicWindow, TE.Us);
TE.sustainedLicks = countEventFromTE(TE, 'Port1In', sustainedWindow, TE.Us);

pupLag = 0.3; % pupil dilation lags cholinergic signal (Nelson and Mooney)
TE.phasicPup = bpCalcPeak_Pupil(TE.pupil, phasicWindow + pupLag, TE.Us);
TE.sustainedPup = bpCalcPeak_Pupil(TE.pupil, sustainedWindow + pupLag, TE.Us);
TE.pupBaseline = bpCalcPeak_Pupil(TE.pupil, [1 4]);

%%
TE.Wheel = processTrialAnalysis_Wheel(sessions, 'duration', 11, 'Fs', 20, 'startField', 'Start');
TE.wheelBaseline = mean(TE.Wheel.data.V(:,bpX2pnt(-3, 20, -4):bpX2pnt(0, 20, -4)), 2);
%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = uigetdir;
% basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\';
% basepath = 'Z:\SummaryAnalyses\LickNoLick_odor_v2_BaselineTrialByTrial\DC_20and17_combined\';
sep = strfind(TE.filename{1}, '_');
subjectName = TE.filename{1}(1:sep(2)-1);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);


%%
truncateSessionsFromTE(TE, 'init', 'usLicks', filterTE(TE, 'trialType', 1));

%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% cross sessions bleaching curve and exponential fits
for channel = channels
    figname = ['sessionBleach_Correction_ch' num2str(channel)];
    ensureFigure(figname, 1);
    plot(TE.Photometry.data(channel).blF_raw, 'k'); hold on;
    plot(TE.Photometry.data(channel).blF, 'r');
    if saveOn
        saveas(gcf, fullfile(savepath, [figname '.fig']));
        saveas(gcf, fullfile(savepath, [figname '.jpg']));
    end
    % cross trial bleaching fits for each session plotted as axis array
    if channel == 1
        figname = ['trialBleach_Correction_ch' num2str(channel)];
        ensureFigure(figname, 1);
        nSessions = size(TE.Photometry.bleachFit, 1);
        subA = ceil(sqrt(nSessions));
        for counter = 1:nSessions
            subplot(subA, subA, counter);
            plot(TE.Photometry.bleachFit(counter, channel).trialTemplateFullX, TE.Photometry.bleachFit(counter, channel).trialTemplateFull, 'g'); hold on;            
            plot(TE.Photometry.bleachFit(counter, channel).trialTemplateX, TE.Photometry.bleachFit(counter, channel).trialTemplate, 'b'); 
            plot(TE.Photometry.bleachFit(counter, channel).fitX, TE.Photometry.bleachFit(counter, channel).trialFit, 'r');
        %     title(num2str(counter));    
        end
        % average of all trials for this channel to eyeball correction
        figname2 = ['corrected_allTrials_ch1']; 
        ensureFigure(figname2, 1);
        [ha, hl] = phPlotAverageFromTE(TE, 1:length(TE.filename), 1,...
    'FluorDataField', 'ZS', 'window', [0.1, max(TE.Photometry.xData) - min(TE.Photometry.xData)], 'zeroTimes', TE.Photometry.startTime); %high value, reward
        if saveOn
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.fig'));
            saveas(gcf, fullfile(savepath, 'trialBleach_Correction.jpg'));
            saveas(gcf, fullfile(savepath, [figname2 '.fig']));
            saveas(gcf, fullfile(savepath, [figname2 '.jpg']));            
        end
    end
end
%% generate trial lookups for different combinations of conditions
LNL_conditions;

%% all behavior CsPlus

    %%
    ensureFigure('all_behavior_CsPlus', 1);
    reversals = find(diff(TE.BlockNumber(csPlusTrials, :))) + 1;
    subplot(1,5,1);
    
    image(TE.Wheel.data.V(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
    set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('Velocity');
    ylabel('trial number');

%     subplot(1,5,2);
%     title('pupil');
%     try
%         image(TE.pupil.pupDiameterNorm(csPlusTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%         set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
%         line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
%         colormap('parula');  
%         title('Pupil Diameter');    
%     catch
%     end
    subplot(1,5,2);
    imagesc(TE.whisk.whisk(csPlusTrials, :), 'XData', [-4 7], [0 8])
    title('whisking');
    
    subplot(1,5,3);
    eventRasterFromTE(TE, csPlusTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    set(gca, 'XLim', [-4 7]);
    title('licking');
    
    
    subplot(1,5,4); phRasterFromTE(TE, csPlusTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    subplot(1,5,5); phRasterFromTE(TE, csPlusTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('DAT');    
axs = findobj(gcf, 'Type', 'axes');

set(axs, 'FontSize', 18);
set(axs(2:end), 'YTick', []);
saveas(gcf, fullfile(savepath, 'allBehavior_whisk_csPlus'), 'fig'); 
saveas(gcf, fullfile(savepath, 'allBehavior_whisk_csPlus'), 'jpeg');
