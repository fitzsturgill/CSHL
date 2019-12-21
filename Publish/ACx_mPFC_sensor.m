load('Z:\SummaryAnalyses\ChAT_PE_Manuscript\ACx_mPFC\ACh_Coherence\ACh_Coherence.mat');
% ACh43_03 is best example
sensor = ACh_Coherence; % easier to type
clear ACh_Coherence;

Fs = str2num(sensor.Param.SamplingRate);

savePath = 'Z:\SummaryAnalyses\ChAT_PE_Manuscript\ACx_mPFC\';
saveOn = 1;
%% construct data matrices for the different animals
include = {'ACh43_03', 'ACh43_05'};
allData = struct(...
    'data', [],...
    'xData', [],...
    'name', cell(1,1),...
    'minPoints', []...
    );
allData = repmat(allData, length(include), 1);
for counter = 1:length(include)
    mouse = include{counter};    
    sd1 = []; sd2 = []; % pain to figure out max number of points
    minPoints = Inf;
    for trial = 1:5
        sd1 = expandVertCat(sd1, sensor.(mouse).(['Trial' num2str(trial)]).Photometry{1}', 'left');
        sd2 = expandVertCat(sd2, sensor.(mouse).(['Trial' num2str(trial)]).Photometry{2}', 'left');
        minPoints = min(minPoints, length(sensor.(mouse).(['Trial' num2str(trial)]).Photometry{1}));
    end
    allData(counter).data = cat(3, sd1, sd2);    
    allData(counter).xData = (0:size(sd1, 2) - 1) / Fs;
    allData(counter).minPoints = minPoints;
    allData(counter).name = mouse;
end

%% example mouse
mouse = 'ACh43_05';
pos = find(strcmp(include, mouse));
minPoints = allData(pos).minPoints;
data1 = allData(pos).data(:,1:minPoints,1)';
data2 = allData(pos).data(:,1:minPoints,2)';

data1_raw = data1;
data2_raw = data2;
% data1 = detrend(data1);
data1 = zscore(data1);
% data2 = detrend(data2);
data2 = zscore(data2);

xData = allData(pos).xData;


%% find good trials


% trial = 5;
figSize = [3 1.5];

trialsToShow = 4;
duration = 120;
Fs = 20    ;
saveName = sprintf('Acx_mPFC_pick_example_traces');
ensureFigure(saveName, 1);
for trial= 1:5

    subplot(5, 1, trial); hold on;

    
    lh = plot(xData(1:size(data1, 1)), [data1(:,trial) data2(:,trial)]);
    lh(1).Color = [1 0 0];
    lh(2).Color = [0 1 0];
end


%% trial 1 18-38s,     other 3, 13-33s


win = [18 38];
trial = 1;
figSize = [3 0.75];

saveName = 'wheel_ACx_mPFC_Figure_exampleTrace';
ensureFigure(saveName, 1);
axes; hold on;

p1 = bpX2pnt(win(1), Fs);
p2 = bpX2pnt(win(2), Fs);
data = zeros(p2 - p1 + 1, 2); %
data(:,1) = data1(p1:p2,trial);
data(:,2) = data2(p1:p2,trial);
data = nanzscore(data);
data = data + [0 0];
xData2 = xData(p1:p2);

% rewX = TE.Reward{trial}(:,1) - TE.PhotometryHF.startTime(trial);
% rewX = rewX((rewX > win(1)) & (rewX < win(end)));
% rewX = rewX - xData(1);
% rewX = repmat(rewX', 2, 1);
xData2 = xData2 - xData2(1);
% rewY = repmat([max(data(:)); min(data(:))], 1, size(rewX, 2));
% plot(rewX, rewY, 'Color', [0 0 1]);

lh = plot(xData2, data, 'k');
lh(1).Color = [0 0 0];
lh(2).Color = [0.5 0.5 0.5];
set(gca, 'YTick', [0 2], 'XTick', [0 10 20], 'XTickLabel', {}, 'YTickLabel', {});

formatFigurePublish('size', figSize);
if saveOn    
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));
end

%% coherence for ACh_43_03



params.Fs = Fs;
params.trialave = 1;
params.err = [2 0.05];
params.tapers = [3 5];
params.pad = 1;
params.fpass = [0.001 10]; % can only go to 10Hz which is the Nyquist frequency

[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data1, data2, params);
% f(1) = eps;


si = randperm(size(data2, 2));
[C_shuff,phi_shuff,~,~,~,f_shuff,~, ~, Cerr_shuff] = coherencyc(data1(:,si), data2, params);


%%
figSize = [1.5 1];

saveName = 'ACx_mPFC_coherence';
ensureFigure(saveName, 1); 
axes; hold on;

f(1) = eps;
subplot(1,1,1); hold on;
hl = [];
[th, ~] = boundedline(f_shuff, C_shuff, Cerr_shuff(1,:)' - C_shuff, 'k');
hl(end + 1) = th;
[th, ~] = boundedline(f, C, Cerr(1,:)' - C, 'b');
hl(end + 1) = th;

set(gca, 'XScale', 'log', 'XLim', [0.01 10], 'YLim', [-0.1 1.25]);
% xlabel('Frequency');
% ylabel('Coherence');
% legend(hl, {'matched', 'shuffled'}, 'Box', 'Off');

formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));
end


%% cross correlation
figSize = [1.5 1];
saveName = 'ACx_mPFC_xCorr';
ensureFigure(saveName, 1);
maxLagInSeconds = 5;
maxLag = round(maxLagInSeconds * Fs);

[r, shiftR, rawR, lags] = correctedXCorr(data1, data2, maxLag);
ensureFigure('xcorr', 1);
plot(lags * 1/Fs, r, 'k'); hold on;
% plot(lags * 1/Fs, rawR, 'r');
% plot(lags * 1/Fs, shiftR, 'b');
% xlabel('Time (s)'); ylabel('mPFC x ACx XCorr'); 
set(gca, 'XTick', [-5 0 5], 'YTick', [0 0.5], 'YLim', [-0.05 0.7]);
% legend({'corrected', 'raw', 'shift predictor'});
        
formatFigurePublish('size', figSize);
if saveOn    
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));
end

%% noise correlation 

% sensor.ACh43_03.RewardFluoFiber1
mouse = 'ACh43_03';
blWin = [1 10];
rWin = [11 30];

phPeak_rew = zeros(size(sensor.(mouse).RewardFluoFiber1, 1), 2);

phPeak_rew(:,1) = mean(sensor.(mouse).RewardFluoFiber1(:, blWin(1):blWin(2)), 2) - mean(sensor.(mouse).RewardFluoFiber1(:, rWin(1):rWin(2)), 2);
phPeak_rew(:,2) = mean(sensor.(mouse).RewardFluoFiber2(:, blWin(1):blWin(2)), 2) - mean(sensor.(mouse).RewardFluoFiber2(:, rWin(1):rWin(2)), 2);

ensureFigure('acx_mpfc_noise_scatter', 1);
scatter(phPeak_rew(:,1), phPeak_rew(:,2));
setXYsameLimit;
% correlation is 0.38
corr(phPeak_rew(:,1), phPeak_rew(:,2))



