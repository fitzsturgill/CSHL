saveOn = 1;
%%
sessions = bpLoadSessions([], 'DC_56_wheel_v1_Sep14_2018_Session3.mat', 'Z:\FitzRig2\Data\DC_56\wheel_v1\Session Data\'); % load sessions
%%
trial = 1;
channel = 1;
Fs = 6100;
rawDataPre = sessions.SessionData.NidaqData{trial,1}(:,channel);
rawDataPre = rawDataPre';
ref = sessions.SessionData.NidaqData{trial,2}; % reference signal parameters

t = 0:1/Fs:(30 - 1/Fs);
refData =  sin(2*pi*ref.freq(channel) .*t + ref.phaseShift(channel));
refData_90 = sin(2*pi*ref.freq(channel) .*t + ref.phaseShift(channel) + pi/2);
 [z,p,k] = butter(5, 25/(Fs/2), 'high');
[sos, g] = zp2sos(z,p,k);
rawData = filtfilt(sos, g, rawDataPre);


%
mixed_0 = rawData .* refData;
mixed_90 = rawData .* refData_90;



%%
[z,p,k] = butter(5, 25/(Fs/2), 'low');
[sos, g] = zp2sos(z,p,k);


filt_0 = filtfilt(sos, g, mixed_0);
filt_90 = filtfilt(sos, g, mixed_90);
demod = (filt_0 .^2 + filt_90 .^2).^(1/2);

%% fft
% masure the amplitude of the reference signal by taking its FFT

% mod_filt or any acquired data specified in seconds should have an even
% number of samples
toDo = {rawData, mixed_0, mixed_90, filt_0, filt_90, demod};


for counter = 1:length(toDo)
    thisData = toDo{counter};
    L = length(thisData);
    n = 2 ^ nextpow2(L);
    Y = fft(thisData, n);
    

    P2 = abs(Y/n);
    P1 = P2(1:n/2 + 1);
    P1(2:end-1) = 2 * P1(2:end-1);

    f = Fs * (0:(n/2))/n;
    if counter == 1
        Pout = zeros(length(toDo), length(P1));
    end
    Pout(counter, :) = P1;
end


%%
ensureFigure('test', 1); axes;
cmap = colormap('hsv');
set(gca, 'Colormap', cmap);
labels = {'rawData', 'mixed0', 'mixed90', 'filt0', 'filt90', 'demod'};
for counter = 1:length(toDo)
    subplot(2,3,counter);
    plot(f, Pout(counter, :));
    set(gca, 'XLim', [0 1000], 'YLim', [0 0.004]);
    set(gca, 'YScale', 'linear');
    title(labels{counter});
end







