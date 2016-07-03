% odor test script  -   measuring odor kinetics and stability with PID

% valve/ trialtype/ outcome 5 = 10% Isoamyl Acetate --
% valve 6 = 10% ethyl tiglate
% valve 7 = 10% cineole


session = bpLoadSession;
session.SessionData = demodulateSession(session.SessionData);
session.analysis = processAnalysis_Photometry(session.SessionData, [], 'zeroField', 'Odor');

%% already ran demodulateSession and processAnalysis_Photomtery
suffix = '_1in10Dilution';
odors = [5:7];
odorNames = {'IAA', 'ET', 'Cin'};
nOdors = length(odors);
Fs = 6100;
startX = -3;
% baseline and peak measurement indices
bl1 = bpX2pnt(-3, Fs, startX);
bl2 = bpX2pnt(0, Fs, startX);

% measure peak from 0.8 --> 1seconds after valve actuation
p1 = bpX2pnt(0.8, Fs, startX);
p2 = bpX2pnt(1, Fs, startX);


% initialize
xData = session.analysis.Photometry.data(odors(1),odors(1)).x(:,1)';
samples = length(xData);
rawData = cell(nOdors,1);
normData = cell(nOdors, 1); % data scaled from baseline to peak as 0 to 1 for averaging to discern onset latency and kinetics of odor pulse
peaks = cell(nOdors, 1); % pid signal at peak
baselines = cell(nOdors, 1);
delta = cell(nOdors, 1); % peak - baseline
normAverages = zeros(nOdors, samples);
normN = zeros(nOdors, 1);
normSEM = zeros(nOdors, samples);

for i = 1:nOdors
    odor = odors(i);
    thisData = session.analysis.Photometry.data(odor,odor).mod(:,:,1);
    rawData{i,1} = thisData;
    peaks{i, 1} = mean(thisData(:,p1:p2), 2);
    baselines{i, 1} = mean(thisData(:,bl1:bl2), 2);
    delta{i, 1} = peaks{i, 1} - baselines{i, 1};
    normData{i, 1} = (thisData - repmat(baselines{i, 1}, 1, samples)) ./ (repmat(peaks{i, 1}, 1, samples) - repmat(baselines{i, 1}, 1, samples));
    normAverages(i, :) = mean(normData{i, 1});
    normN(i, :) = size(thisData, 1); % number of trials
    normSEM(i, :) = std(normData{i, 1}) / sqrt(normN(i,:));
end
    
    
% peaks figure
ensureFigure('Peaks', 1); hold on
for i=1:nOdors
    plot(delta{i,1});
end
legend('IAA', 'ET', 'Cin');
saveas(gcf, ['Peaks' suffix '.fig']);
saveas(gcf, ['Peaks' suffix '.jpg']);

% averages figure
ensureFigure('Averages', 1); hold on
plot(repmat(xData', 1, nOdors), normAverages');
legend('IAA', 'ET', 'Cin');
saveas(gcf, ['Averages' suffix '.fig']);
saveas(gcf, ['Averages' suffix '.jpg']);

% allData figure
ensureFigure('rawData', 1); 
for i=1:nOdors
    subplot(ceil(sqrt(nOdors)),ceil(sqrt(nOdors)),i); % make it square with room to spare :)
    plot(rawData{i,1}');
    title(odorNames{i});
end
saveas(gcf, ['rawData' suffix '.fig']);
saveas(gcf, ['rawData' suffix '.jpg']);
    


clear('ans');
save(['odorTestScript' suffix]);

