% wheel_exampleMouse


DB = dbLoadExperiment('wheel');
saveOn = 1;


% load the session without reward

% load('Z:\SummaryAnalyses\wheel_v1\DC_56_wheel_v1_Aug27_2018_Session1\TE.mat');
animal = 'DC_56';
dbLoadAnimal(DB, animal);
savepath = fullfile(DB.path, [filesep 'figure' filesep]);
ensureDirectory(savepath);

%% examples

% saveName = '

%% coherence


Photometry = 'PhotometryHF';


data_chat = TE.(Photometry).data(1).ZS';
% data_chat = diff(data_chat, 1, 1);
% data_chat = data_chat(:,2:end) - data_chat(:,1:end-1); % whiten
% data_chat = nanzscore2(data_chat); % standardize
% data_chat = nanzscore(data_chat); % standardize
data_dat = TE.(Photometry).data(2).ZS';
% data_dat = diff(data_dat, 1, 1);
% data_dat = data_dat(:,2:end) - data_dat(:,1:end-1); % whiten
% data_dat = nanzscore2(data_dat); % standardize
% data_dat = nanzscore(data_dat); % standardize


params.Fs = TE.(Photometry).sampleRate;
duration = length(TE.(Photometry).xData) / params.Fs;
params.trialave = 1;
params.err = [2 0.05];
params.tapers = [5 9];
params.pad = 1;
params.fpass = [1/duration 20];

[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat, data_dat, params);
Cerr = max(Cerr, 0);
f(1) = eps; % for log scale, you can't show zero



% [C_scram,phi_scram,S12_scram,S1_scram,S2_scram,f_scram,confC_scram, phistd_scram, Cerr_scram] = coherencyc(data_chat(:,si), data_dat, params);
% f_scram(1) = eps;

[C_scram,phi_scram,S12_scram,S1_scram,S2_scram,f_scram,confC_scram, phistd_scram, Cerr_scram] = coherencyc(circshift(data_chat, 1, 2), data_dat, params);
Cerr_scram = max(Cerr_scram, 0);
f_scram(1) = eps;

%%
fpass_avg = [0.2 3]; 
figSize = [1.5 1];
saveName = 'wheel_coherence';
ensureFigure(saveName, 1);
% x = linspace(1,100, length(f_scram));
% y = randn(length(f_scram), 1) + 10;
% b = ones(size(Cerr_scram)) .* 4;
axes; hold on;
% boundedline(x, y, b', 'k');%, 'alpha');
% abs(cohData(counter).cohS.Cerr' - cohData(counter).cohS.C)
% boundedline(f, C_scram, abs(Cerr_scram' - C_scram), 'r');%, 'alpha');
% boundedline(f, C, abs(Cerr' - C), 'k');%, 'alpha');


% pdf not written properly with 2 bounded line series for some odd reason,
% so just using bounded line for the matched, not trial-shifted coherence
plot(repmat(fpass_avg, 2, 1), repmat([0; 1], 1, 2), 'k--');
plot(f_scram, C_scram, 'k');
plot(f_scram, Cerr_scram', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.25);
boundedline(f, C, abs(Cerr' - C), 'b');%, 'alpha');
% plot(f_scram, C_scram, 'k');%, 'alpha');
% plot(f_scram, Cerr_scram', 'Color', [0.8 0.8 0.8]);
% plot(f, C, 'b');%, 'alpha');
% plot(f, Cerr', 'Color', [0 0 0.8]);

set(gca, 'XScale', 'log', 'XLim', params.fpass);
set(gca, 'XTick', [0.1 1 10], 'YLim', [0 1]);

xlabel('Frequency');
ylabel('Coherence');
%
formatFigurePublish('size', figSize);
if saveOn    
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
end


%% example traces

% #5 is the best, shows nice reward responses, coherence, and
% arousal/movement modulation
% good ones potentially: 2,5,8,10,11


trial = 5;
figSize = [3 1.5];


saveName = sprintf('wheel_example_traces_FINAL_tr%d', trial);

ensureFigure(saveName, 1);

ax = axes; hold on;
Fs = 20;
duration = 120;

data = zeros(Fs*duration, 5); % wheel pupil whisk chat dat
data(:,1) = nanzscore(TE.Wheel.data.V(trial,1:Fs*duration));
data(:,2) = nanzscore(TE.Whisk.whiskNorm(trial,1:Fs*duration));
data(:,3) = nanzscore(TE.pupil.pupDiameterNorm(trial,1:Fs*duration));
data(:,4) = nanzscore(TE.Photometry.data(1).ZS(trial,:));
data(:,5) = nanzscore(TE.Photometry.data(2).ZS(trial,:));
xdata = TE.Photometry.xData;

% splay the data vertically
splay = range(data);
splay = cumsum(splay);
data = data - splay;

rewX = TE.Reward{trial}(:,1) - TE.Photometry.startTime(trial);
rewX = repmat(rewX', 2, 1);
rewY = repmat([max(data(:)); min(data(:))], 1, size(rewX, 2));

plot(rewX, rewY, 'Color', [0.7 0.7 0.7]);
lh = plot(xdata, data, 'k');
lh(4).Color = mycolors('chat');
lh(5).Color = mycolors('dat');
ax.YAxis.Visible = 'off';
ax.XLim = [0 120];
ax.YLim = [min(data(:)) max(data(:))];
xlabel('time (s)');

formatFigurePublish('size', figSize);

if saveOn    
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
end


%%
% 1) Coherence w.r.t. reward for each animal, bar graphs

loadpath = fullfile(DB.path, filesep, 'pooled', filesep);
load(fullfile(loadpath, 'coherence_pooled.mat'), 'cohData');
disp(['*** loaded: ' fullfile(loadpath, 'coherence_pooled.mat') ' ***']);



%% now the bar graph:
% PARAMETERS MUST MATCH THOSE IN WHEEL_POOLED_COHERENCE.M
figSize = [1.7 1];
mix = find(strcmp(animal, DB.animals)); % mouse index
fpass_avg = [0.2 3]; 
bins = [0 5; 5 120]; % 0-5s = peri,   5-120s = post
movingWin = [5 0.5]; % 10 0.5
maxTime = 40;
saveName = 'wheel_coherence_ex_timeFromReward';
ensureFigure(saveName, 1);

% bin means
if fpass_avg(1)
    ix1 = crossing(cohData(mix).cxcg.f - fpass_avg(1));
else
    ix1 = 1; % in case you are starting from frequency 0
end
ix2 = crossing(cohData(mix).cxcg.f - fpass_avg(2));
bins = [0:movingWin(2):maxTime];
[xC_Means, xC_Errors, xC_timeFromReward] = binnedMeansXY(cohData(mix).tfr, mean(cohData(mix).C(ix1:ix2,:))', bins);
[xC_Means_shift, xC_Errors_shift, ~] = binnedMeansXY(cohData(mix).tfr, mean(cohData(mix).CS(ix1:ix2,:))', bins);

axes; hold on;
% errorbar(xC_timeFromReward, xC_Means, xC_Errors, 'b');
% errorbar(xC_timeFromReward, xC_Means_shift, xC_Errors_shift, 'k');
boundedline(xC_timeFromReward, xC_Means_shift, xC_Errors_shift, 'k');
boundedline(xC_timeFromReward, xC_Means, xC_Errors, 'b');


xlabel('time from reward(s)'); ylabel('Coh. 0.2-3Hz');


formatFigurePublish('size', figSize);

if saveOn    
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));
end

                 
%%




















% 
% animal = 'DC_56';
% 
% 
% dbLoadAnimal(DB, animal);
% 
% 
% % cross correlogram from pre-reward period
% 
% 
% s1 = extractDataByTimeStamps(TE.Photometry.data(1).ZS, TE.Photometry.startTime, 20, TE.Reward, [-20 0]);
% s2 = extractDataByTimeStamps(TE.Photometry.data(2).ZS, TE.Photometry.startTime, 20, TE.Reward, [-20 0]);


%% example traces, all in a pdf:

% #5 is the best, shows nice reward responses, coherence, and
% arousal/movement modulation
% good ones potentially: 2,5,8,10,11


trials = 1:length(TE.filename);
fh = zeros(length(trials), 1);
for counter = 1:length(trials)
    trial = trials(counter);
    saveName = sprintf('wheel_example_traces_tr%d', trial);
    fh(counter) = ensureFigure(saveName, 1);
    mcLandscapeFigSetup(gcf);
    ax = axes; hold on;
    Fs = 20;
    duration = 120;

    data = zeros(Fs*duration, 5); % wheel pupil whisk chat dat
    data(:,1) = nanzscore(TE.Wheel.data.V(trial,1:Fs*duration));
    data(:,2) = nanzscore(TE.Whisk.whiskNorm(trial,1:Fs*duration));
    data(:,3) = nanzscore(TE.pupil.pupDiameterNorm(trial,1:Fs*duration));
    data(:,4) = nanzscore(TE.Photometry.data(1).ZS(trial,:));
    data(:,5) = nanzscore(TE.Photometry.data(2).ZS(trial,:));
    xdata = TE.Photometry.xData;

    % splay the data vertically
    splay = range(data);
    splay = cumsum(splay);
    data = data - splay;

    rewX = TE.Reward{trial}(:,1) - TE.Photometry.startTime(trial);
    rewX = repmat(rewX', 2, 1);
    rewY = repmat([max(data(:)); min(data(:))], 1, size(rewX, 2));

    plot(rewX, rewY, 'Color', [0.7 0.7 0.7]);
    lh = plot(xdata, data, 'k');
    lh(4).Color = mycolors('chat');
    lh(5).Color = mycolors('dat');
    ax.YAxis.Visible = 'off';
    ax.XLim = [0 120];
    ax.YLim = [min(data(:)) max(data(:))];
    xlabel('time (s)');
    textBox(sprintf('Trial #%d', trial));
end


% write pdfs
h = waitbar(0, 'slowly writing pdfs');
pdfn = fullfile(savepath, 'wheel_example_traces_all.pdf');
for counter = 1:length(fh)
        if counter == 1
        export_fig(fh(counter),pdfn);  % write to pdf
    else
        export_fig(fh(counter),'-append',pdfn);  % write to pdf
    end
    waitbar(counter/length(fh));
end
close(h);


