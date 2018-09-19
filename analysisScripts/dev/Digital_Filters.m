    %% try filtering first
    % note-   5 pole Butterworth filter in Matlab used in Frohlich and McCormick  
     % Create butterworth filter
    lowCutoff = 15/6100 * 2; % multiply by 2 to convert to rad/sample- see butter documentation
    % for a cutoff freq of 300Hz and sample rate of 1000Hz, cuto2e12ff
    % corresponds to 0.6pi rad/sample    300/1000 * 2 = 0.6    
    [b, a] = butter(10, lowCutoff, 'low'); % mimic filt filt in which filter order is doubled...
    
    [z,p,k] = butter(10,2 * pi * 15, 'low', 's');
    [sos, g] = zp2sos(z,p,k);
    [b, a] = zp2tf(z,p,k);
    [h, w] = freqs(b, a, 2^12);
    ensureFigure('filters', 1);  
    plot(w/(15 * 2 * pi), mag2db(abs(h)));
%     hfvt = fvtool(b,a,sos);
%     legend(hfvt, 'TF Design', 'ZPK Design');
    
    %%
    [b2, a2] = cheby1(10,3,lowCutof(f));
    
    
    
    
%%     last thing on desktop today


%% try filtering first
% note-   5 pole Butterworth filter in Matlab used in Frohlich and McCormick
% Create butterworth filter
%    lowCutoff = 15/6100 * 2; % multiply by 2 to convert to rad/sample normalized to the nyquest frequency
%    - see butter documentation
% for a cutoff freq of 300Hz and sample rate of 1000Hz, cuto2e12ff
% corresponds to 0.6pi rad/sample    300/1000 * 2 = 0.6


[z,p,k] = butter(10,15/(6100/2), 'low'); % normalize to Nyquist frequency, gives units of pi radians per sample in terms of the Nyquist frequency
[b2, a2] = butter(10,15/(6100/2), 'low'); % normalize to Nyquist frequency, gives units of pi radians per sample in terms of the Nyquist frequency
%    [sos, g] = zp2sos(z,p,k);
[b, a] = zp2tf(z,p,k);
w = linspace(0.1/(6100/2), 50/(6100/2), 1000);
f = linspace(0.1, 50, 1000);
[h] = freqs(b, a, w);
[h2] = freqs(b2, a2, w);
ensureFigure('filters', 1);
plot(f, mag2db(abs(h)), 'b'); hold on;
plot(f, mag2db(abs(h2)), 'r--')

%%  try without making w (freq vector) yourself, now using 's' form where you specify cutoff frequency in radians
[z,p,k] = butter(30,15 * 2 * pi, 'low', 's'); 
%    [sos, g] = zp2sos(z,p,k);
[b, a] = zp2tf(z,p,k);
[b2, a2] = butter(10, 15 * 2 * pi, 'low', 's');
% [b,a] = butter(10,15/(6100/2), 'low');
% w = linspace(0.1/(6100/2), 50/(6100/2), 1000);
% f = linspace(0.1, 50, 1000);
[h, w] = freqs(b, a, 2^12);
[h2, w2] = freqs(b2, a2, 2^12);
ensureFigure('filters', 1);
% fitz = w .* (6100/2/w(end));
plot(w/(2 * pi), mag2db(abs(h)), 'b'); hold on;
plot(w2/(2 * pi), mag2db(abs(h2)), 'r--');
%% are 


lowCutoff = 15/(6100/2); % normalize to Nyquist frequency, gives units of pi radians per sample in terms of the Nyquist frequency

[b, a] = butter(10, lowCutoff, 'low'); % mimic filt filt in which filter order is doubled...

[z,p,k] = butter(10, lowCutoff, 'low');
[sos, g] = zp2sos(z,p,k);


ensureFigure('filters_tool', 1);

    hfvt = fvtool(b, a, sos, 'FrequencyScale', 'log');
    legend(hfvt, 'TF Design', 'ZPK Design');
    
    
    
    
%%
n = 6;
Wn = [2.5e6 29e6]/500e6;
ftype = 'bandpass';

% Transfer Function design
[b,a] = butter(n,Wn,ftype);      % This is an unstable filter

% Zero-Pole-Gain design
[z,p,k] = butter(n,Wn,ftype);
sos = zp2sos(z,p,k);

% Display and compare results
hfvt = fvtool(b,a,sos,'FrequencyScale','log');
legend(hfvt,'TF Design','ZPK Design')
    
    