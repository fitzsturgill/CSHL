
% Fitz Sturgill 2018
%% make dummy data 
 

Fs_raw = 6100;
Fs = 610;


duration = 10;
t = linspace(0, duration, duration * Fs_raw);
f1 = 4; % 5hz for ch1 and 2
f2 = 10;
fmod1 = 211;
fmod2 = 531;
refAmp = 0.6;
sigAmp = 0.1;
sig1 = sin(2*pi*f1 * t) * sigAmp + (sigAmp + 0.1); % signals with offset to mimic real calcium data
sig2 = sin(2*pi*f2 * t) * sigAmp + (sigAmp + 0.1);   % signals with offset to mimic real calcium data

ref1 = sin(2*pi*fmod1 * t);
ref2 = sin(2*pi*fmod2 * t);

mod = sig1 .* (ref1 * refAmp + refAmp) + sig2 .* (ref2 * refAmp + refAmp);
    mod = mod + 0.5; % add an offset to mimic dark current/voltage of photodetector

figure; subplot(3,1,1);   plot(t, [sig1' + 0.1, sig2' + 0.1]); set(gca, 'XLim', [0 2], 'YLim', [0 2]); legend('signal1', 'signal2'); title('raw signals');
subplot(3,1,2); plot(t, [ref1' * refAmp + refAmp, ref2' * refAmp + refAmp]); set(gca, 'XLim', [0 0.1]); title('reference signals'); legend('ref1', 'ref2');
subplot(3,1,3); plot(t, mod, 'g');  set(gca, 'XLim', [0 2]); set(gca, 'YLim', [0 2]); title('modulated data'); xlabel('time (s)');
 
 
 %% AC filter modulated data and plot filtered version of signal
 figure;
 [z,p,k] = butter(5, 25/(Fs_raw/2), 'high'); % 25Hz cutoff frequency normalized to Nyquist frequency
[sos, g] = zp2sos(z,p,k);
mod_filt = filtfilt(sos, g, mod);
plot(t, mod_filt , 'g');  set(gca, 'XLim', [0 2]);
title('AC filtered data');

%% demodulate the synthetic data and compare to original prior to modulation
figure;
subplot(2,1,1); title('signal 1');
% make lowpass filter
[z,p,k] = butter(5, 15/(Fs_raw/2), 'low'); % 15Hz cutoff frequency normalized to Nyquist frequency
[sos, g] = zp2sos(z,p,k);

ref_0 = ref1;
ref_90 = sin(2*pi*fmod1 * t + pi/2); % phase shift by 90 degrees
mixed_0 = mod_filt .* ref_0;
mixed_90 = mod_filt .* ref_90;

mixed_0_filt = filtfilt(sos, g, mixed_0);
mixed_90_filt = filtfilt(sos, g, mixed_90);

demod = (mixed_0_filt .^2 + mixed_90_filt .^2) .^(1/2); % take pythagorean distance
demod = demod * 2 / refAmp;
% demod = demod * 2;

 plot(t, sig1, 'b'); hold on; plot(t, demod, 'g--'); set(gca, 'XLim', [0 2], 'YLim', [0 2]);
legend('original', 'demodulated');

% signal 2
subplot(2,1,2); title('signal 2');
% ref_0 = ref2;
ref_0 = sin(2*pi*fmod2 * t); % phase shift by 90 degrees
ref_90 = sin(2*pi*fmod2 * t + pi/2); % phase shift by 90 degrees
mixed_0 = mod_filt .* ref_0;
mixed_90 = mod_filt .* ref_90;

mixed_0_filt = filtfilt(sos, g, mixed_0);
mixed_90_filt = filtfilt(sos, g, mixed_90);

demod = (mixed_0_filt .^2 + mixed_90_filt .^2) .^(1/2); % take pythagorean distance
demod = demod * 2 / refAmp;
% demod = demod * 2;

 plot(t, sig2, 'r'); hold on; plot(t, demod, 'y--'); set(gca, 'XLim', [0 2], 'YLim', [0 2]);
legend('original', 'demodulated');
%% FFT of modulated data



Y = fft(mod_filt);
L = length(mod_filt);
P2 = abs(Y/length(mod_filt));
P1 = P2(1:L/2 + 1);
P1(2:end-1) = 2 * P1(2:end-1);
f = Fs_raw * (0:(L/2))/L;

ensureFigure('fft', 1); plot(f, P1);
set(gca, 'XLim', [0 1000]); xlabel('frequency (Hz)');



