% 10/15/17
% research_statement_PID_vs_Phasic

pidPath = 'Z:\SummaryAnalyses\PID\odorTest\Session Data';
pidFile = 'Dummy Subject_odorTest_Mar28_2016_Session3.mat'; % file for 1/10 odor dilutions


% decimate and load PID data into a matrix, generate vector of trial types, and xData 
sessions = bpLoadSessions([], pidFile, pidPath);
%%
decimationFactor = 50;
Fs = 6100 / decimationFactor;
nTrials = sessions.SessionData.nTrials;
nSamples = size(sessions.SessionData.NidaqData{1,1}(:,3), 1); % PID is in 3rd channel
pidData = zeros(nTrials, nSamples/decimationFactor);
trialTypes = sessions.SessionData.TrialTypes(1:nTrials);
trialTypes = trialTypes(:);
for counter = 1:nTrials
    pidData(counter, :) = decimate(sessions.SessionData.NidaqData{counter,1}(:,3)', decimationFactor);
end
% trial type key:  5 = Amyl Acetate, 6 = Ethyl Tiglate, 7 = Cineole 

xData = 0:1/Fs:6 - 1/Fs;
xData = xData - 3; % zero to odor onset time

pidAvgs = zeros(3, size(pidData, 2));
pidAvgs(1, :) = mean(pidData(trialTypes == 5, :));
pidAvgs(2, :) = mean(pidData(trialTypes == 6, :));
pidAvgs(3, :) = mean(pidData(trialTypes == 7, :));
% normalize
pidAvgs = bsxfun(@minus, pidAvgs, mean(pidAvgs(:, 1:bpX2pnt(0, Fs, -3)), 2));
pidAvgs = bsxfun(@rdivide, pidAvgs, max(pidAvgs, [], 2));
% slope
pidSlope = diff(pidAvgs, 1, 2);
pidSlope = bsxfun(@rdivide, pidSlope, max(pidSlope, [], 2));
ensureFigure('test', 1); plot(xData, pidAvgs', '-*'); hold on;
plot(xData(2:end), pidSlope', '--*');
set(gca, 'XLim', [-0.2 0.5])
%% run grand averages script to load CuedOutcome_Odor_Complete summary averages
cuedOutcome_pooledAnalysis_grandAverages;

%%
% note! In Uchida Mainen, 2003,   they measure the field response in the
% olfactary epithelium.   The onset (measured as the 2nd derivative of
% voltage) was 172ms, they subtracted off a component to yield a more
% conservative (or lenient depending on your persepective as a speed or
% integration person) estimate of the lag as 140ms.
% I think for PID measurements 1st derivative makes sense since you first
% derivative of voltage yields current and second derivative yields current
% slope

% lets start by assuming a lag of 200ms from eyeballing the derivative of
% the PID traces
ensureFigure('phasic_speed');
pm = [2, 4];
for counter = 1:7 
    subplot(pm(1), pm(2), counter);
    % first chunk (3rd dimension) is high value
    plot(xData_ph - 0.2, cuePh_norm(counter, :, 1));
    set(gca, 'XLim', [-0.2 0.5]);
end



window = [-0.5 0.5];
pidStart = -3;
phStart = -3;

% 






