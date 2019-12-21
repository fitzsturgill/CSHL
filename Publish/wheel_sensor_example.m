
saveOn = 1;
savePath = 'Z:\SummaryAnalyses\ChAT_PE_Manuscript\wheel_sensor_example\';
load('Z:\SummaryAnalyses\wheel_v1\ACh_3_wheel_v1_Feb04_2019_Session1\TE.mat');



%% find good trials


trial = 5;
figSize = [3 1.5];

trialsToShow = 4;
duration = 120;
Fs = TE.PhotometryHF.sampleRate;
for counter = 1:3
    trials = (1:trialsToShow) + (counter-1)*trialsToShow;
    saveName = sprintf('wheel_ex_tr%dto%d', trials(1), trials(end));
    ensureFigure(saveName, 1);
    for tcounter = 1:trialsToShow
        trial = trials(tcounter);
        if trial > length(TE.filename)
            break
        end
        subplot(trialsToShow, 1, tcounter); hold on;
        
        data = zeros(Fs*duration, 2); %
        data(:,1) = nanzscore(TE.PhotometryHF.data(1).ZS(trial,:));
        data(:,2) = nanzscore(TE.PhotometryHF.data(2).ZS(trial,:));
        xdata = TE.PhotometryHF.xData;
        rewX = TE.Reward{trial}(:,1) - TE.PhotometryHF.startTime(trial);
        rewX = repmat(rewX', 2, 1);
        rewY = repmat([max(data(:)); min(data(:))], 1, size(rewX, 2));
        plot(rewX, rewY, 'Color', [0 0 0]);

        lh = plot(xdata, data);
        lh(1).Color = [1 0 0];
        lh(2).Color = [0 1 0];
    end
end

%% use trial 5, 40s - 80s into trial

win = [60 80];
trial = 5;
figSize = [3 0.75];

saveName = 'wheel_sensor_Figure_exampleTrace';
ensureFigure(saveName, 1);
axes; hold on;

p1 = bpX2pnt(win(1), Fs);
p2 = bpX2pnt(win(2), Fs);
data = zeros(p2 - p1 + 1, 2); %
data(:,1) = TE.PhotometryHF.data(1).ZS(trial,p1:p2);
data(:,2) = TE.PhotometryHF.data(2).ZS(trial,p1:p2);
data = nanzscore(data);
data = data + [0 0];
xData = TE.PhotometryHF.xData(p1:p2);

rewX = TE.Reward{trial}(:,1) - TE.PhotometryHF.startTime(trial);
rewX = rewX((rewX > win(1)) & (rewX < win(end)));
rewX = rewX - xData(1);
rewX = repmat(rewX', 2, 1);
xData = xData - xData(1);
rewY = repmat([max(data(:)); min(data(:))], 1, size(rewX, 2));
plot(rewX, rewY, 'Color', [0 0 1]);

lh = plot(xData, data, 'k');
lh(1).Color = [0 0 0];
lh(2).Color = [0.5 0.5 0.5];
set(gca, 'YTick', [0 2], 'XTick', [0 10 20], 'XTickLabel', {}, 'YTickLabel', {});

formatFigurePublish('size', figSize);
if saveOn    
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));
end

%% 

Photometry = 'PhotometryHF';
data1 = TE.(Photometry).data(1).ZS';
data2 = TE.(Photometry).data(2).ZS';
% data1 = nanzscore(data1);
% data2 = nanzscore(data2);



params.Fs = Fs;
params.trialave = 1;
params.err = [2 0.05];
params.tapers = [3 5];
params.pad = 1;
params.fpass = [0.001 20];

[C,~,~,~,~,f,~,~,Cerr] = coherencyc(data1, data2, params);
% f(1) = eps;

% scramble trial labels
si = randperm(size(data2, 2));
[C_shuff,phi_shuff,~,~,~,f_shuff,~, ~, Cerr_shuff] = coherencyc(data1(:,si), data2, params);
% f_shuff(1) = eps;

%% coherence graph
figSize = [1.5 1];

saveName = 'wheel_sensor_example_BLA_coherence';
ensureFigure(saveName, 1); 
axes; hold on;

f(1) = eps;
subplot(1,1,1); hold on;
hl = [];
[th, ~] = boundedline(f, C, Cerr(1,:)' - C, 'b');
hl(end + 1) = th;
[th, ~] = boundedline(f_shuff, C_shuff, Cerr_shuff(1,:)' - C_shuff, 'k');
hl(end + 1) = th;
set(gca, 'XScale', 'log', 'XLim', [0.01 20], 'YLim', [-0.1 1.25]);
% xlabel('Frequency');
% ylabel('Coherence');
% legend(hl, {'matched', 'shuffled'}, 'Box', 'Off');

formatFigurePublish('size', figSize);
if saveOn 
    print(gcf, '-dpdf', fullfile(savePath, [saveName '.pdf']));
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));
end

        