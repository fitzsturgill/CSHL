sessions = bpLoadSessions([], 'ACh_6_wheel_v1_Feb23_2019_Session1.mat', 'Z:\FitzRig2\Data\ACh_6\wheel_v1\Session Data\');
load('Z:\SummaryAnalyses\wheel_405test\ACh_6_wheel_v1_Feb23_2019_Session1\TE.mat')

%%
data_470 = TE.Photometry_ch1.data(1).raw';
data_470 = data_470(:);

data_405 = TE.Photometry_ch2.data(1).raw';
data_405 = data_405(:);

nPoints = size(TE.Photometry_ch1.data(1).raw, 2);
nTrials = size(TE.Photometry_ch1.data(1).raw, 1);
xData = repmat([NaN 1:nPoints-1] ./ TE.Photometry_ch1.sampleRate, nTrials, 1);


trialStarts = sessions.SessionData.TrialStartTimestamp' - sessions.SessionData.TrialStartTimestamp(1);

xData = xData + trialStarts;
xData = xData';
xData = xData(:);

%%
ensureFigure('weird_405effect', 1);
plot(xData, data_470, 'b'); hold on;
plot(xData, data_405, 'm'); 
xlabel('time (s)');
legend('470nm decoded', '405nm decoded');