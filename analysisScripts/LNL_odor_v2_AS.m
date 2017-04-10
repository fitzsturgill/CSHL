
saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_LNL_odor_V2(sessions);
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial');


%% extract peak trial dFF responses to cues and reinforcement and lick counts
csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
TE.phPeakMean_cs = bpCalcPeak_dFF(TE.Photometry, 1, csWindow, TE.Cue, 'method', 'mean');
TE.csLicks = countEventFromTE(TE, 'Port1In', csLickWindow, TE.Us);

usWindow = [0 0.5];
usZeros = cellfun(@(x,y,z,a) max(x(1), max(y(1), max(z(1), a(1)))), TE.Reward, TE.Punish, TE.WNoise, TE.Neutral); %'Reward', 'Punish', 'WNoise', 'Neutral'
TE.phPeakMean_us = bpCalcPeak_dFF(TE.Photometry, 1, [0 0.5], usZeros, 'method', 'mean');
TE.phPeakPercentile_us = bpCalcPeak_dFF(TE.Photometry, 1, [0 0.5], usZeros, 'method', 'percentile', 0.9);



TE.usLicks = countEventFromTE(TE, 'Port1In', [0 2], TE.Us);

csLicks = 
%%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
basepath = uigetdir;
subjectName = TE.filename{1}(1:7);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);