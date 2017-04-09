
saveOn = 1;
%% 
sessions = bpLoadSessions;
%%
TE = makeTE_CuedOutcome_Odor_Complete(sessions);
TE.Photometry = processTrialAnalyszis_Photometry2(sessions, 'dFFMode', 'expFit', 'blMode', 'byTrial');


%% extract peak trial dFF responses to cues and reinforcement and lick counts
csWindow = cellfun(@(x,y,z) [x(1) max(y(end), z(end))], TE.Cue, TE.AnswerLick, TE.AnswerNoLick); % max- to select either AnswerLick or AnswerNoLick timestamp (unused state contains NaN)
TE.phPeak_cs = bpCalcPeak_dFF(TE.Photometry, 1, csWindow, TE.Cue, 'method', 'mean');
TE.phPeak_us = bpCalcPeak_dFF(TE.Photometry, 1, [0 0.5], TE.Us, 'method', 'mean');


TE.csLicks = countEventFromTE(TE, 'Port1In', csLickWindow, TE.Us);
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