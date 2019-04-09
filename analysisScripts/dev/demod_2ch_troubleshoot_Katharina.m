



sessions = bpLoadSessions([], 'T01_DetectionConfidence_Mar29_2019_Session2.mat', 'Z:\FitzRig1\Data\Katharina\');

%% make TE
nTrials = sessions.SessionData.nTrials;
Fs = sessions.SessionData.NidaqData{1,2}.sample_rate(1);
TE = struct(...
    'LED1_amp', zeros(nTrials, 1),...
    'LED2_amp', zeros(nTrials, 1),...
    'duration', zeros(nTrials, 1),...
    'nSamples', zeros(nTrials, 1)...
    );
for counter = 1:nTrials
    TE.LED1_amp(counter) = sessions.SessionData.NidaqData{counter,2}.amp(1);
    TE.LED2_amp(counter) = sessions.SessionData.NidaqData{counter,2}.amp(1);
    TE.duration(counter) = size(sessions.SessionData.NidaqData{counter,1}, 1) / Fs;
    TE.duration(counter) = size(sessions.SessionData.NidaqData{counter,1}, 1);
end

%% demodulate
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'simple', 'blMode', 'byTrial',...
    'channels', [1 2], 'baseline', [0 3], 'downsample', 305, 'forceAmp', 1, 'uniformOutput', 0, 'startField', 'StartRecording');

%%

ensureFigure('demodulated', 1);
ylim = [0 .2];
subplot(3,2,1);
plot(TE.Photometry.data(1).raw{7}); set(gca, 'YLim', ylim);
ylabel('blue');
subplot(3,2,2);
plot(TE.Photometry.data(2).raw{7}); set(gca, 'YLim', ylim);
subplot(3,2,3);
plot(TE.Photometry.data(1).raw{8}); set(gca, 'YLim', ylim);
subplot(3,2,4);
plot(TE.Photometry.data(2).raw{8}); set(gca, 'YLim', ylim);
subplot(3,2,5);
plot(TE.Photometry.data(1).raw{9}); set(gca, 'YLim', ylim);
subplot(3,2,6);
plot(TE.Photometry.data(2).raw{9}); set(gca, 'YLim', ylim);

%%
ensureFigure('unmodulated', 1);
ylim = [0 .2];
subplot(3,2,1);
plot(TE.Photometry.data(1).raw{4}); set(gca, 'YLim', ylim);
subplot(3,2,2);
plot(TE.Photometry.data(2).raw{4}); set(gca, 'YLim', ylim);
subplot(3,2,3);
plot(TE.Photometry.data(1).raw{5}); set(gca, 'YLim', ylim);
subplot(3,2,4);
plot(TE.Photometry.data(2).raw{5}); set(gca, 'YLim', ylim);
subplot(3,2,5);
plot(TE.Photometry.data(1).raw{6}); set(gca, 'YLim', ylim);
subplot(3,2,6);
plot(TE.Photometry.data(2).raw{6}); set(gca, 'YLim', ylim);
