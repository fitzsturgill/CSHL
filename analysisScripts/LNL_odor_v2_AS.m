% Not going to change processTrialAnalysis_Photometry2 just yet...  But I need to change how I deal with channels, 
% and especially zeroField- to make it flexible for different timings (as in operant conditioning) 
saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_LNL_odor_V2(sessions);
% assume that photometry channels are consistent across sessions
channels=[];
if sessions.SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
end

if sessions.SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
end

TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial', 'zeroField', 'Cue', 'channels', channels);


%% extract peak trial dFF responses to cues and reinforcement and lick counts
% csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
nTrials = length(TE.filename);
csWindow = zeros(nTrials, 2);
csWindow(:,1) = cellfun(@(x) x(1), TE.Cue);
csWindow(:,2) = cellfun(@(x,y) max(x(end), y(end)), TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
csWindow = bsxfun(@minus, csWindow, TE.Photometry.startTime);



TE.csLicks = countEventFromTE(TE, 'Port1In', csWindow, TE.Cue);

usWindow = [0 0.5];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'

for channel = channels
    TE.phPeakMean_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'mean');
    TE.phPeakPercentile_cs(channel) = bpCalcPeak_dFF(TE.Photometry, channel, csWindow, TE.Cue, 'method', 'percentile', 'percentile', 0.9);
    TE.phPeakMean_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.5], usZeros, 'method', 'mean');
    TE.phPeakPercentile_us(channel) = bpCalcPeak_dFF(TE.Photometry, channel, [0 0.5], usZeros, 'method', 'percentile', 'percentile', 0.9);
end



TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], usZeros);

%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = uigetdir;
subjectName = TE.filename{1}(1:7);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);


%%
truncateSessionsFromTE(TE, 'init');
%%
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end

%% cross sessions bleaching curve and dual exponential fits
for channel = channels
    figname = ['sessionBleach_Correction_ch' num2str(channel)];
    ensureFigure(figname, 1);
    plot(TE.Photometry.data(channel).blF_raw, 'k'); hold on;
    plot(TE.Photometry.data(channel).blF, 'r');
    if saveOn
        saveas(gcf, fullfile(savepath, [figname '.fig']));
        saveas(gcf, fullfile(savepath, [figname '.jpg']));
    end
    %% cross trial bleaching fits for each session plotted as axis array
    figname = ['trialBleach_Correction_ch' num2str(channel)];
    ensureFigure(figname, 1);
    nSessions = length(TE.Photometry.bleachFit);
    subA = ceil(sqrt(nSessions));
    for counter = 1:nSessions
        subplot(subA, subA, counter);
        plot(TE.Photometry.bleachFit(counter).trialTemplate, 'k'); hold on;
        plot(TE.Photometry.bleachFit(counter).trialFit, 'r');
    %     title(num2str(counter));    
    end
    if saveOn
        saveas(gcf, fullfile(savepath, 'trialBleach_Correction.fig'));
        saveas(gcf, fullfile(savepath, 'trialBleach_Correction.jpg'));
    end
end

