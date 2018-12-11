saveOn = 0;
%%
saveOn = 1;
%% load eiother sessions 11,12 or 13
sessions = bpLoadSessions([], 'test_offset_wheel_v1_Nov30_2018_Session11.mat', 'Z:\FitzRig1\Data\test_offset\wheel_v1\Session Data\'); % load sessions
%% 
TE = makeTE_wheel_v1(sessions); % make TE structure
channels = [1];
dFFMode{1} = 'simple';
baselineEnd = sessions(1).SessionData.NidaqData{1,2}.duration(1) - 0.2;
        
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', 'simple', 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', 1, 'baseline', [1 baselineEnd], 'startField', 'Baseline', 'downsample', 10, 'ACfilter', 1);
Fs_raw = 6100;

%%

t = linspace(0, 6, 6 * Fs_raw);
ensureFigure('compare_mod_and_demod', 1);
ylims = [-5 10];
subplot(2,2,1);
plot(TE.Photometry.xData(610:end-610), TE.Photometry.data(1).raw(3,610:end-610));
subplot(2,2,2);
plot(TE.Photometry.xData(610:end-610), TE.Photometry.data(1).raw(4,610:end-610));

% plot(t(6100:end), sessions.SessionData.NidaqData{3,1}(6100:end,1));

subplot(2,2,3);
plot(t(6100:6710), sessions.SessionData.NidaqData{3,1}(6100:6710,1));
subplot(2,2,4);
plot(t(6100:6710), sessions.SessionData.NidaqData{4,1}(6100:6710,1));
% set(gca, 'YLim', ylims);

%%

%% spectral profile of demodulated and raw data
    regularTrials = rem(TE.trialNumber, 2) == 0;
    offsetTrials = ~regularTrials;    
    Fs = TE.Photometry.sampleRate;    
    freqRange = [0 50];    
    nSamples = length(TE.Photometry.xData);    
    deltaX = 1/Fs_raw;       
     
    [z,p,k] = butter(5, 25/(6100/2), 'high');
    [sos, g] = zp2sos(z,p,k);
    
%     regularRawData = zeros(sum(regularTrials), 36600);
%     offsetRawData = zeros(sum(offsetTrials), 36600);
    allData = zeros(sessions.SessionData.nTrials, 36600);
    for counter = 1:sessions.SessionData.nTrials
        allData(counter, :) = filtfilt(sos, g, sessions.SessionData.NidaqData{counter,1})';
%             allData(counter, :) = sessions.SessionData.NidaqData{counter,1}';
    end
    regularRawData = allData(regularTrials, :);
    offsetRawData = allData(offsetTrials, :);
    % lets just plot the power spectra
    params.fpass = freqRange;
    params.Fs = Fs;
    params.trialave = 1;
    powerFig = ensureFigure('powerFigure', 1);
    subplot(2,2,1, 'YScale', 'linear', 'XScale', 'linear'); hold on        
    [S, f] = mtspectrumc(TE.Photometry.data(1).raw(regularTrials, :)', params);
    plot(f, S , 'g'); 
    [S, f] = mtspectrumc(TE.Photometry.data(1).raw(offsetTrials, :)', params);
    plot(f, S, 'k');    
    title('GCaMP channel, demodulated'); ylabel('Power');
    legend('regular', 'linear range');
    
    params.Fs = 6100;
    params.fpass = [0 1000];
    subplot(2,2,2, 'YScale', 'log', 'XScale', 'linear'); hold on        

    [S, f] = mtspectrumc(offsetRawData', params);
    plot(f, S, 'k-');    
    [S, f] = mtspectrumc(regularRawData', params);
    plot(f, S , 'g-'); 
    title('GCaMP channel, raw but AC filtered'); ylabel('Power');
    legend('linear range', 'regular');

    subplot(2,2,3, 'YScale', 'log', 'XScale', 'linear'); hold on        
    [S, f] = mtspectrumc(offsetRawData', params);
    plot(f, S, 'k-');    
    title('GCaMP channel, raw but AC filtered');
    subplot(2,2,4, 'YScale', 'log', 'XScale', 'linear'); hold on        
    [S, f] = mtspectrumc(regularRawData', params);
    plot(f, S , 'g-'); 
    title('GCaMP channel, raw but AC filtered'); ylabel('Power');