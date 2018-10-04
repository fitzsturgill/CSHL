%% Examples showing how to use phAlignedWindow and phAverageFromTE (new version)


%% example 1
% constant zeros/ windows
trial_subset = csPlusTrials; % this can be any subset of trials
trial_subset2 = csMinusTrials;
[alignedPhotometryData, xData] = phAlignedWindow(TE, trial_subset, 1, 'window', [-3 5], 'zeroTimes', TE.Cue);
ensureFigure('constantWin', 1);
subplot(1,2,1);
imagesc('XData', xData', 'CData', alignedPhotometryData);
subplot(1,2,2);
avgData = phAverageFromTE(TE, trial_subset, 1, 'window', [-3 5], 'zeroTimes', TE.Cue);
plot(avgData.xData, avgData.Avg);

%% example 2
%  you can also supply zero times as a vector rather than a cell array
[alignedPhotometryData, xData] = phAlignedWindow(TE, trial_subset, 1, 'window', [-4 5], 'zeroTimes', TE.Photometry.startTime + 3.456789); % zeroed 3.456789 seconds after photometry start
ensureFigure('arbitraryZero', 1); imagesc('XData', xData', 'CData', alignedPhotometryData);


%% example 3
% here AnswerLick would only be defined for 'hit' trials so the output data
% would have lines of just NaNs if the animal didn't lick during the answer
% period
myzeros = cellfun(@(x) x(1), TE.PostUsRecording);
winStart = cellfun(@(x) x(1), TE.AnswerLick) - myzeros;
winEnd = zeros(size(myzeros)) + 2; % 2 seconds after
mywin = [winStart winEnd];
[alignedPhotometryData, xData] = phAlignedWindow(TE, trial_subset, 1, 'window', mywin, 'zeroTimes', myzeros); 
ensureFigure('variableWindow1', 1); 
imagesc('XData', xData', 'CData', alignedPhotometryData);






%% example 4
% test variable zero/alignment points
TE.firstLick = calcEventLatency(TE, 'Port1In', TE.AnswerStart, TE.Us); % my calcEventLatency functions computes the latency between Bpod time stamps and events (e.g. licks and the start of a bpod state)
TE.firstLick = TE.firstLick + cellfun(@(x) x(1), TE.AnswerStart); % convert to bpod time (relative to bpod trial start)

% test a variable window and varaible zeros
varWin = [-4 + rand(length(TE.filename), 1) 3 + rand(length(TE.filename), 1)];
[alignedPhotometryData, xData] = phAlignedWindow(TE, trial_subset, 1, 'window', varWin, 'zeroTimes', TE.firstLick);
ensureFigure('variableWindow2', 1); 
subplot(1,3,1);
imagesc('XData', xData, 'CData', alignedPhotometryData);
avgData = phAverageFromTE(TE, {trial_subset, csMinusTrials}, 1, 'window', varWin, 'zeroTimes', TE.firstLick);

subplot(1,3,2); plot(avgData.xData', avgData.Avg'); title('averages');
subplot(1,3,3); plot(avgData.xData', avgData.N'); ylabel('number of components in average at each time point');

