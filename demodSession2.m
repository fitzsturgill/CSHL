
lowCutoff = 30; % optional
ch = 1;
trial = 102;
sampleRate = 6100; % req
modRate = 211; % req, get out of trial settings don't hard codee

rawData = session.SessionData.NidaqData{trial,1}(:,ch); %req
nSamples = length(rawData);
refData = session.SessionData.NidaqData{trial,2}(1:nSamples,ch);
tBaseline = session.SessionData.TrialSettings(trial).PreTrialRecording; % optional? or req?

preSamples = floor(tBaseline / (1/sampleRate));


[processedData, ~, ~] = zscoreByRange(rawData, 1, preSamples);
[refData, ~, ~] = zscore(refData);

% generate 90degree shfited copy of refData
samplesPerPeriod = 1/modRate / (1/sampleRate);
quarterPeriod = round(samplesPerPeriod / 4); % ideally you shouldn't have to round, i.e. mod frequencies should be close to factors of sample freq
refData90 = circshift(refData, [1 quarterPeriod]);

processedData_0 = processedData .* refData;
processedData_90 = processedData .* refData90;

demodData = (processedData_0 .^2 + processedData_90 .^2) .^(1/2); % quadrature decoding

 % Create butterworth filter
lowCutoff = lowCutoff/sampleRate; % originally (mcFilter) I multiplied by a factor of 2, why did I do this?

% note-   5 pole Butterworth filter in Matlab used in Frohlich and McCormick  
[b, a] = butter(5, lowCutoff, 'low');         
demodDataFilt = filtfilt(b,a,demodData);

fig = ensureFigure('test', 1);
plot(demodDataFilt);

