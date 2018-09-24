saveOn = 0;
%%
saveOn = 1;
%%
sessions = bpLoadSessions([], 'DC_56_wheel_v1_Sep14_2018_Session3.mat', 'Z:\FitzRig2\Data\DC_56\wheel_v1\Session Data\'); % load sessions
%% 
TE = makeTE_wheel_v1(sessions); % make TE structure

channels=[]; dFFMode = {}; %BL = {};
if sessions(1).SessionData.Settings.GUI.LED1_amp > 0
    channels(end+1) = 1;
    dFFMode{end+1} = 'simple';
%     BL{end + 1} = [1 4];
end

if sessions(1).SessionData.Settings.GUI.LED2_amp > 0
    channels(end+1) = 2;
    dFFMode{end+1} = 'simple';
%     BL{end + 1} = [2 4];    
end

baselineEnd = sessions(1).SessionData.NidaqData{1,2}.duration(1) - 0.2;
        
TE.Photometry = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [1 baselineEnd], 'startField', 'Baseline', 'downsample', 10, 'ACfilter', 1);

TE.PhotometryDC = processTrialAnalysis_Photometry2(sessions, 'dFFMode', dFFMode, 'blMode', 'byTrial',...
    'zeroField', 'Baseline', 'channels', channels, 'baseline', [1 baselineEnd], 'startField', 'Baseline', 'downsample', 10, 'ACfilter', 0);
Fs_raw = 6100;

%%
 [z,p,k] = butter(5, 25/(Fs_raw/2), 'low');
[sos, g] = zp2sos(z,p,k);
t = linspace(0,30, 30 * Fs_raw);
ensureFigure('compare_mod', 1);
ylims = [-5 10];
subplot(2,2,1);
plot(TE.Photometry.xData(610:end-610), TE.Photometry.data(1).raw(3,610:end-610));
subplot(2,2,2);
plot(t(6100:end), sessions.SessionData.NidaqData{3,1}(6100:end,1));

subplot(2,2,3);
plot(TE.PhotometryDC.xData(610:end), TE.PhotometryDC.data(1).raw(4,610:end)); 
subplot(2,2,4);
plot(t(6100:end), filtfilt(sos, g, sessions.SessionData.NidaqData{4,1}(6100:end,1)));
% set(gca, 'YLim', ylims);

%% eyeball the amplitudes of the unmodulated and modulated data
Fs_raw = 6100;
Fs = 610;
% make dummy data and demodulate it to check why the amplitudes don't match between DC and AC modes....

% I think the answer is that 1) you need to subtract off the dark
% current/voltage from the photodetector  2) (check how big these are given
% that I saw them on log/log plot but) you are probably not recovering the
% side lobes that seem to appear on power spectra of modulated data around
% the reference frequency
t = linspace(0, 30, 30 * Fs_raw);
f1 = 5; % 5hz for ch1 and 2
fmod = 531;
refAmp = 0.2;
testData = sin(2*pi*f1 * t) * 0.06 + .4;  % 0.02 amplitude signal with 0.3 amplitude offset
testRef = sin(2*pi*fmod * t) * refAmp + refAmp;  % reference 
mod = testData .* testRef + 0.3; % add offset here to make it look like actual data

ensureFigure('eyeball', 1);  plot(t, sessions.SessionData.NidaqData{1,1}(:,1)); hold on; plot(t, mod, 'g');  plot(t, sessions.SessionData.NidaqData{2,1}(:,1), 'r'); legend('modulated', 'unmodulated');
 legend('modulated', 'synthetic', 'unmodulated');
 
 %% AC filter and plot versions of the modulated and test data
 ensureFigure('AC_test', 1); 
 [z,p,k] = butter(5, 25/(Fs_raw/2), 'high');
[sos, g] = zp2sos(z,p,k);
plot(t, filtfilt(sos, g, sessions.SessionData.NidaqData{1,1}(:,1)), 'b'); hold on; 
plot(t, filtfilt(sos, g, mod), 'g'); legend('mod_AC', 'synth_AC');

%% demodulate the synthetic data and compare to original prior to modulation
 [acz,acp,ack] = butter(5, 25/(Fs_raw/2), 'high');
[sos_ac, g_ac] = zp2sos(acz,acp,ack);
mod_filt = filtfilt(sos_ac, g_ac, mod);
% hack
% testData = sin(2 * pi * 5 * t) + 5;
% mod = sin(2 * pi * 5 * t) .* sin(2*pi*fmod*t);
ref_0 = sin(2*pi*fmod * t);
ref_90 = sin(2*pi*fmod * t + pi/2);
mixed_0 = mod_filt .* ref_0;
mixed_90 = mod_filt .* ref_90;
% make filter
[z,p,k] = butter(5, 15/(Fs_raw/2), 'low');
[sos, g] = zp2sos(z,p,k);
mixed_0_filt = filtfilt(sos, g, mixed_0);
mixed_90_filt = filtfilt(sos, g, mixed_90);
testData_filtered = filtfilt(sos, g, testData);
demod = (mixed_0_filt .^2 + mixed_90_filt .^2) .^(1/2); % take pythagorean distance
demod = demod * 2 / refAmp;
% demod = mixed_0_filt * 2 / 0.3;

ensureFigure('synth_demodulation', 1); plot(t, testData_filtered, 'b'); hold on; plot(t, demod, 'r--');
legend('original', 'demodulated');

%% masure the amplitude of the reference signal by taking its FFT

% mod_filt or any acquired data specified in seconds should have an even
% number of samples

Y = fft(mod_filt);
L = length(mod_filt);
P2 = abs(Y/length(mod_filt));
P1 = P2(1:L/2 + 1);
P1(2:end-1) = 2 * P1(2:end-1);
f = Fs_raw * (0:(L/2))/L;
ensureFigure('fft', 1); plot(f, P1);


%% signal to noise
acTrials = rem(TE.trialNumber, 2) > 0;
dcTrials = ~acTrials;

window = [-2 2];
blSamples = (0 - window(1)) * Fs;
clims = [-5 5];

[rewards_dat_ac, ts_ac, tn_ac] = extractDataByTimeStamps(TE.Photometry.data(2).dF(acTrials, :), TE.Photometry.startTime(acTrials), 610, TE.Reward(acTrials), [-2 2]);
rewards_chat_ac = extractDataByTimeStamps(TE.Photometry.data(1).dF(acTrials, :), TE.Photometry.startTime(acTrials), 610, TE.Reward(acTrials), [-2 2]);
rewards_ac = cat(3, rewards_chat_ac, rewards_dat_ac);
[rewards_dat_dc, ts_dc, tn_dc] = extractDataByTimeStamps(TE.PhotometryDC.data(2).dF(dcTrials, :), TE.Photometry.startTime(dcTrials), 610, TE.Reward(dcTrials), [-2 2]);
rewards_chat_dc = extractDataByTimeStamps(TE.PhotometryDC.data(1).dF(dcTrials, :), TE.Photometry.startTime(dcTrials), 610, TE.Reward(dcTrials), [-2 2]);
rewards_dc = cat(3, rewards_chat_dc, rewards_dat_dc);


ensureFigure('phRasters', 1); subplot(2,2,1); imagesc([-2 2], [1 size(rewards_chat_ac, 1)],rewards_chat_ac, clims); set(gca, 'XTick', []); title('ChAT'); ylabel('freq modulated');
subplot(2,2,2); imagesc([-2 2], [1 size(rewards_dat_ac, 1)], rewards_dat_ac, clims); set(gca, 'XTick', []); title('DAT');
subplot(2,2,3); imagesc([-2 2], [1 size(rewards_chat_dc, 1)], rewards_chat_dc, clims); xlabel('time from reward (s)'); ylabel('unmodulated');
subplot(2,2,4); imagesc([-2 2], [1 size(rewards_dat_dc, 1)],rewards_dat_dc, clims); 

%% calculate discriminability of reward responses from the baseline

bl_ac = squeeze(mean(rewards_ac(:, bpX2pnt(-0.5, Fs, -2):bpX2pnt(0, Fs, -2),:)));
peak_ac = squeeze(mean(rewards_ac(:,bpX2pnt(0, Fs, -2):bpX2pnt(0.5, Fs, -2),:)));

bl_dc = squeeze(mean(rewards_dc(:, bpX2pnt(-0.5, Fs, -2):bpX2pnt(0, Fs, -2),:)));
peak_dc = squeeze(mean(rewards_dc(:,bpX2pnt(0, Fs, -2):bpX2pnt(0.5, Fs, -2),:)));

D = []; P = []; CI = [];
for counter = 1:2    
    [D(1,counter), P(1, counter), CI(1,counter,1:2)] = rocarea_CI(peak_ac(:,counter), bl_ac(:,counter), 'boot', 10000, 'scale');
    [D(2,counter), P(2, counter), CI(2,counter,1:2)] = rocarea_CI(peak_dc(:,counter), bl_dc(:,counter), 'boot', 10000, 'scale');
end


ensureFigure('dPrime', 1); 
errors = abs(CI - D);
errorbar([1 1; 2 2], D, squeeze(errors(:,:,1)), squeeze(errors(:,:,1)));
set(gca, 'XLim', [0.5 2.5], 'XTick', [1 2], 'XTickLabel', {'AC', 'DC'});
set(gca, 'YLim', [0 1]); ylabel('Dprime');
legend('ChAT', 'DAT');

%% just plot some unnormalized data
ensureFigure('rawData', 1);
subplot(1,2,1);
plot(TE.Photometry.xData(Fs:end-Fs), TE.Photometry.data(1).raw(3,Fs:end-Fs), 'b'); hold on; 
plot(TE.PhotometryDC.xData(Fs:end-Fs), TE.PhotometryDC.data(1).raw(4,Fs:end-Fs), 'r'); 
% rawX = linspace(1, 1 + 0.033, 201);
rawX = linspace(1, 4, 201);
plot(rawX, sessions.SessionData.NidaqData{3,1}(6100:6300, 1), 'b');
legend('AC', 'DC'); set(gca, 'YLim', [0 5]);
subplot(1,2,2);
plot(TE.Photometry.xData(Fs:end-Fs), TE.Photometry.data(2).raw(3,Fs:end-Fs), 'b'); hold on;
plot(TE.PhotometryDC.xData(Fs:end-Fs), TE.PhotometryDC.data(2).raw(4,Fs:end-Fs), 'r'); 
plot(rawX, sessions.SessionData.NidaqData{3,1}(6100:6300, 2), 'b');
% plot(repmat(min
legend('AC', 'DC');set(gca, 'YLim', [0 5]);
%% spectral profile of demodulated data
    
    
    Fs = TE.Photometry.sampleRate;
    
    freqRange = [0 100];
    
    nSamples = length(TE.Photometry.xData);
    
    deltaX = 1/Fs_raw;
    
    
    colors = {'y', 'r', 'k', 'b', 'g', 'm'};

    % lets just plot the spectrograms of the first acqs

    params.fpass = freqRange;
    params.Fs = Fs;
    params.trialave = 1;
    powerFig = ensureFigure('powerFigure', 1);
    axes('YScale', 'log', 'XScale', 'linear'); hold on
        
    [S, f] = mtspectrumc(TE.Photometry.data(1).ZS(acTrials, :)', params);
    plot(f, S , 'g'); 
    [S, f] = mtspectrumc(TE.Photometry.data(2).ZS(acTrials, :)', params);
    plot(f, S, 'r');
    [S, f] = mtspectrumc(TE.PhotometryDC.data(1).ZS(dcTrials, :)', params);
    plot(f, S, 'y');    
    [S, f] = mtspectrumc(TE.PhotometryDC.data(2).ZS(dcTrials, :)', params);
    plot(f, S, 'm');    
    

    xlabel('Frequency');
    ylabel('Power');
    title('Noise Spectra');
        
    
