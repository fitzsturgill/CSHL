
trial = 1;
% t = linspace(0, 10, 1000);
h = ensureFigure('lockin', 1);

% plot(t, cos(t) + 2, 'k'); hold on;
% plot(t, (cos(t) + 2) .* (cos(t) + 2), 'r');

rawData = session.SessionData.NidaqData{trial,1}(:,1);
nSamples = size(rawData, 1);
refData = session.SessionData.NidaqData{trial,2}(1:nSamples,1);
sampleRate = 6100;  % need to not hardcode in acq software

modF = 211;
[finalData, dmf, dmp] = decode_lockin_fn(rawData, refData, [], modF, sampleRate, 0);

xData = linspace(0, nSamples / sampleRate, nSamples);

plot(xData, finalData);

