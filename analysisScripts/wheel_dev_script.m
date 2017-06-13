%%
 rawTimes = sessions(1).SessionData.RawEvents.Trial{1}.Events.Port3In;
 wheelTimes = rawTimes(rawTimes < 30);
 wheelY = ones(size(wheelTimes));
 wheelTimes = [0 wheelTimes 30];
 wheelY = cumsum([0 wheelY 0]);
 
 [wheelY_new, wheelTimes_new] = resample(wheelY, wheelTimes, 20, 'linear');
 
 ensureFigure('test', 1);
 plot(wheelTimes, wheelY); hold on;
 plot(wheelTimes_new, wheelY_new);
 
 %%
D = 12.7; % diameter of wheel in cm
pr = 200; % pulses/rotation 
dpp = pi * 12.7 / 200;

wheelTimes = rawTimes(rawTimes < 30);
edges = 0:1/20:30;
% wheel_X = cumsum(histcounts(rawTimes, edges) * dpp);
wheel_X = smooth(cumsum(histcounts(rawTimes, edges)), 20);  
wheel_t = 0:1/20:30-1/20;
  
ensureFigure('test', 1);
subplot(2,1,1);
plot(wheelTimes, cumsum(ones(size(wheelTimes)) * dpp), 'b'); hold on;
plot(wheel_t, wheel_X, 'g');

subplot(2,1,2);
plot(wheel_t(1:end-1), diff(wheel_X), 'b'); hold on
plot(wheel_t, gradient(wheel_X), 'g');

%%
window = [-2 2];
[rewards_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2]);
rewards_chat = extractDataByTimeStamps(TE.Photometry.data(1).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2]);

% local dFF
bl_dat = nanmean(rewards_dat(:,1:40), 2);
rewards_dat = bsxfun(@minus, rewards_dat, bl_dat);
rewards_dat = bsxfun(@rdivide, rewards_dat, bl_dat);
bl_chat = nanmean(rewards_chat(:,1:40), 2);
rewards_chat = bsxfun(@minus, rewards_chat, bl_chat);
rewards_chat = bsxfun(@rdivide, rewards_chat, bl_chat);

nTrials = length(sessions.SessionData.nTrials);
ts_abs = zeros(size(ts));
for counter = 1:length(ts)
    ts_abs(counter) = ts(counter) + sessions.SessionData.TrialStartTimestamp(tn(counter));    
end

iri_pre = [Inf; diff(ts_abs)];
iri_post = [diff(ts_abs); Inf];

[~, I] = sort(iri_pre);
iri_pre_sorted = iri_pre(I);
rewards_dat_sorted = rewards_dat(I, :);
rewards_chat_sorted = rewards_chat(I, :);

ensureFigure('random_rewards', 1); 
% subplot(3,2,1); imshow(rewards_chat_sorted, [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); colormap('jet'); title('chat'); ylabel('sorted');
% subplot(3,2,2); imshow(rewards_dat_sorted, [min(min(rewards_dat_sorted)), max(max(rewards_dat_sorted))]); colormap('jet'); title('dat');
subplot(3,2,1); image(rewards_chat, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); colormap('jet');  title('ChAT');
subplot(3,2,2); image(rewards_dat, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', [min(min(rewards_dat_sorted)), max(max(rewards_dat_sorted))]); colormap('jet');    title('DAT');
xdata = linspace(window(1), window(2), size(rewards_chat, 2));
subplot(3,2,3); plot(xdata, nanmean(rewards_chat));
subplot(3,2,4); plot(xdata, nanmean(rewards_dat));
subplot(3,2,5); triggeredEventRasterFromTE(TE, 'Port1In', TE.Reward);

if saveOn
    saveas(gcf, fullfile(savepath, 'random_rewards.fig'));
    saveas(gcf, fullfile(savepath, 'random_rewards.jpg'));
end

%% coherence

data_chat = TE.Photometry.data(1).raw';
% data_chat = data_chat(:,2:end) - data_chat(:,1:end-1); % whiten
% data_chat = nanzscore2(data_chat); % standardize
data_chat = nanzscore(data_chat); % standardize
data_dat = TE.Photometry.data(2).raw';
% data_dat = data_dat(:,2:end) - data_dat(:,1:end-1); % whiten
% data_dat = nanzscore2(data_dat); % standardize
data_dat = nanzscore(data_dat); % standardize


params.Fs = 20;
params.trialave = 1;
params.err = [2 0.05];
params.tapers = [3 5];
params.pad = 1;

[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat, data_dat, params);
f(1) = eps;
ensureFigure('coherence', 1);
% boundedline(f, C, Cerr
subplot(1,1,1); %plot(f,C, 'r'); hold on;
boundedline(f, C, Cerr(1,:)' - C, 'b', 'alpha')
% plot(f, Cerr(1,:), 'm');
% plot(f, Cerr(2,:), 'm');

%
% subplot(1,2,2); plot(f, phi, 'r'); hold on;
% plot(f, phi + 2 * phistd, 'm');
% plot(f, phi - 2 * phistd, 'm');
%
% scramble trial labels
si = randperm(size(data_dat, 2));

[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat(:,si), data_dat, params);
f(1) = eps;
subplot(1,1,1); %plot(f,C, 'k'); hold on;
boundedline(f, C, Cerr(1,:)' - C, 'alpha', 'k');
set(gca, 'XScale', 'log', 'XLim', [0.01 10]);
xlabel('Frequency');
ylabel('Coherence');
formatFigureGRC;
if saveOn
    saveas(gcf, fullfile(savepath, 'coherence.fig'));
    saveas(gcf, fullfile(savepath, 'coherence.jpg'));
    saveas(gcf, fullfile(savepath, 'coherence.epsc'));
end




%%

[rewards_wheel, ts, tn] = extractDataByTimeStamps(TE.Wheel.data.V, TE.Wheel.startTime, 20, TE.Reward, [-1 1]);
ensureFigure('rewards_wheel', 1);
plot(nanmean(rewards_wheel));

%% plot individual trials with reward times annotated

channel = 2;
ensureFigure(['annotated_' num2str(channel)], 1);
for trial = 1:18;

    subplot(6,3,trial);
    trial = trial + 18 * 1;
    ydata = TE.Photometry.data(channel).raw(trial, :);    
    plot(TE.Photometry.xData, ydata, 'k'); hold on;
    tsx = repmat(TE.Reward{trial}(:,1), 1, 2)';
    tsy = [repmat(min(min(ydata)), 1, size(tsx, 2)); repmat(max(max(ydata)), 1, size(tsx, 2))];    
    plot(tsx, tsy, 'r');
end

%% cross correlation
trial = 1;
maxLagInSeconds = 10;
Fs = 20;
maxLag = round(maxLagInSeconds * Fs);
% [r, lags] = xcorr(data_chat(:,trial), data_dat(:,trial), maxLag);
[r, shiftR, rawR, lags] = correctedXCorr(data_chat, data_dat, maxLag);
ensureFigure('xcorr', 1);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('ChAT x DAT XCorr'); 
legend({'corrected', 'raw', 'shift predictor'});


if saveOn
    saveas(gcf, fullfile(savepath, 'xcorr.fig'));
    saveas(gcf, fullfile(savepath, 'xcorr.jpg'));
end

%% pupil data
%  [wheelY_new, wheelTimes_new] = resample(wheelY, wheelTimes, 20, 'linear');

TE = addPupilometryToTE(TE, 'duration', 30, 'zeroField', 'Baseline', 'startField', 'Baseline', 'frameRate', 60, 'frameRateNew', 20);

%%

maxLagInSeconds = 10;
Fs = 20;
maxLag = round(maxLagInSeconds * Fs);
% [r, lags] = xcorr(data_chat(:,trial), data_dat(:,trial), maxLag);

h = ensureFigure('xcorr_pupil', 1);
mcPortraitFigSetup(h);

subplot(3,2,1);
[r, shiftR, rawR, lags] = correctedXCorr(data_chat, data_dat, maxLag);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('ChAT x DAT XCorr'); 
% legend({'corrected', 'raw', 'shift predictor'}, 'Location', 'northwest');

subplot(3,2,2);
[r, shiftR, rawR, lags] = correctedXCorr(data_chat, data_chat, maxLag);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('ChAT AutoCorr'); 
legend({'corrected', 'raw', 'shift predictor', 'Location', 'northwest'});

subplot(3,2,3);
[r, shiftR, rawR, lags] = correctedXCorr(data_chat, data_chat, maxLag);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('DAT AutoCorr'); 
% legend({'corrected', 'raw', 'shift predictor'});

data_pupil = TE.pupil.pupDiameter';
% data_pupil = data_pupil(:,2:end) - data_pupil(:,1:end-1); % whiten
% data_pupil = nanzscore2(data_pupil);
data_pupil = nanzscore(data_pupil);
maxLagInSeconds = 10;
Fs = 20;
maxLag = round(maxLagInSeconds * Fs);
% [r, lags] = xcorr(data_chat(:,trial), data_dat(:,trial), maxLag);
[r, shiftR, rawR, lags] = correctedXCorr(data_chat, data_pupil, maxLag);
subplot(3,2,4);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('ChAT x pupil XCorr'); 
% legend({'corrected', 'raw', 'shift predictor'});

[r, shiftR, rawR, lags] = correctedXCorr(data_dat, data_pupil, maxLag);
subplot(3,2,5);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('DAT x pupil XCorr'); 
% legend({'corrected', 'raw', 'shift predictor'});

[r, shiftR, rawR, lags] = correctedXCorr(data_pupil, data_pupil, maxLag);
subplot(3,2,6);
plot(lags * (1/20), r, 'k'); hold on;
plot(lags * (1/20), rawR, 'r');
plot(lags * (1/20), shiftR, 'b');
xlabel('Time (s)'); ylabel('Pupil AutoCorr'); 
% legend({'corrected', 'raw', 'shift predictor'});



if saveOn
    saveas(gcf, fullfile(savepath, 'xcorr_pupil.fig'));
    saveas(gcf, fullfile(savepath, 'xcorr_pupil.jpg'));
end
%%
%% coherence with pupil

validTrials = find(sum(isnan(data_pupil)) == 0);

params.Fs = 20;
params.trialave = 1;
params.err = [2 0.1];
params.tapers = [3 5];


[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat(:,validTrials), data_pupil(:,validTrials), params);
ensureFigure('coherence_pupil', 1);
subplot(1,2,1); plot(f,C, 'r'); hold on;
plot(f, Cerr(1,:), 'm');
plot(f, Cerr(2,:), 'm');
%
% scramble trial labels
si = randperm(length(validTrials));
si2 = validTrials(si);

[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat(:,si2), data_pupil(:,validTrials), params);

plot(f,C, 'b'); 
plot(f, Cerr(1,:), 'c');
plot(f, Cerr(2,:), 'c');
set(gca, 'XScale', 'log');
xlabel('Frequency');
ylabel('Coherence, chat vs pupil');


[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_dat(:,validTrials), data_pupil(:,validTrials), params);
subplot(1,2,2); plot(f,C, 'r'); hold on;
plot(f, Cerr(1,:), 'm');
plot(f, Cerr(2,:), 'm');

% scramble trial labels
si = randperm(size(data_dat, 2));

[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_dat(:,si2), data_pupil(:,validTrials), params);

plot(f,C, 'b'); 
plot(f, Cerr(1,:), 'c');
plot(f, Cerr(2,:), 'c');
set(gca, 'XScale', 'log');
xlabel('Frequency');
ylabel('Coherence, dat vs pupil');

if saveOn
    saveas(gcf, fullfile(savepath, 'coherence_pupil.fig'));
    saveas(gcf, fullfile(savepath, 'coherence_pupil.jpg'));
end

%% scatter
ensureFigure('scatter_pupil', 1); 
subplot(2,2,1); scatter(reshape(data_chat, numel(data_chat), 1), reshape(data_pupil, numel(data_pupil), 1), '.');
ylabel('pupil'); xlabel('chat');
subplot(2,2,2); scatter(reshape(data_dat, numel(data_dat), 1), reshape(data_pupil, numel(data_pupil), 1), '.');
ylabel('pupil'); xlabel('dat');
subplot(2,2,3); scatter(reshape(TE.Wheel.data.V, numel(data_dat), 1), reshape(data_pupil, numel(data_pupil), 1), '.');
ylabel('pupil'); xlabel('velocity');
%% coherence bar graph


data_chat = TE.Photometry.data(1).raw';
% data_chat = data_chat(:,2:end) - data_chat(:,1:end-1); % whiten
% data_chat = nanzscore2(data_chat); % standardize
data_chat = nanzscore(data_chat); % standardize
data_dat = TE.Photometry.data(2).raw';
% data_dat = data_dat(:,2:end) - data_dat(:,1:end-1); % whiten
% data_dat = nanzscore2(data_dat); % standardize
data_dat = nanzscore(data_dat); % standardize


params.Fs = 20;
params.trialave = 1;
params.err = [2 0.05];
params.tapers = [3 5];
params.pad = 1;

tf = 0.1; % test frequency
cohAll_C = zeros(3,2); % chat x dat, chat x pupil, dat x pupil
cohAll_confC = zeros(3,2);

% chat vs dat
[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat, data_dat, params);
tfp = nearest(f, tf);
cohAll_C(1,1) = C(tfp);
cohAll_confC(1,1) = confC;
si = randperm(size(data_dat, 2));
[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat(:,si), data_dat, params);
cohAll_C(1,2) = C(tfp);
cohAll_confC(1,2) = confC;

% chat vs pupil
validTrials = find(sum(isnan(data_pupil)) == 0);
[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat(:,validTrials), data_pupil(:,validTrials), params);
cohAll_C(2,1) = C(tfp);
cohAll_confC(2,1) = confC;
si = randperm(length(validTrials));
si2 = validTrials(si);
[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat(:,si2), data_pupil(:,validTrials), params);
cohAll_C(2,2) = C(tfp);
cohAll_confC(2,2) = confC;
% dat vs pupil
[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_dat(:,validTrials), data_pupil(:,validTrials), params);
cohAll_C(3,1) = C(tfp);
cohAll_confC(3,1) = confC;
si = randperm(size(data_dat, 2));
[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_dat(:,si2), data_pupil(:,validTrials), params);
cohAll_C(3,2) = C(tfp);
cohAll_confC(3,2) = confC;


%%
% figure; axes; hold on;
xe = 1:3;
offset = 0.1;
ensureFigure('coherence_barGraph', 1); hold on;
bar([1 3 5], cohAll_C(:,1), .3,  'b');
bar([2 4 6], cohAll_C(:,2), .3, 'k');
legend({'matched trials', 'scrambled trials'});
errorbar([1 3 5], cohAll_C(:,1), cohAll_confC(:,1), 'b.', 'linewidth', 2);
errorbar([2 4 6], cohAll_C(:,2), cohAll_confC(:,2), 'k.', 'linewidth', 2);

set(gca, 'TickDir', 'out', 'YLim', [0 1]);
formatFigureGRC;
ylabel('Coherence (0.1Hz)');
if saveOn
    saveas(gcf, fullfile(savepath, 'coherence_barGraph.fig'));
    saveas(gcf, fullfile(savepath, 'coherence_barGraph.jpg'));
    saveas(gcf, fullfile(savepath, 'coherence_barGraph.epsc'));
end

%% 
xdata = TE.Wheel.xData;
mult = 6;
trials = 1:5;
trials = trials + 5 * mult;
data_wheel = TE.Wheel.data.V';
figName = ['pickTrials_' num2str(trials(1)) ' +'];
ensureFigure('pickTrials', 1);
phax = 1:3:13;
puax = phax + 1;
vax = phax + 2;
for counter = 1:5
    trial = trials(counter);
    subplot(5, 3, phax(counter)); 
    plot(xdata, data_chat(:,trial), 'g'); hold on;
    plot(xdata, data_dat(:,trial), 'r'); set(gca, 'XLim', [0 30]); ylabel(['trial ' num2str(trial)]); 
    if counter == 1
        title('photometry');
    end
    subplot(5, 3, puax(counter));
    plot(xdata, data_pupil(:, trial));set(gca, 'XLim', [0 30]); 
    if counter == 1
        title('pupil');
    end
    subplot(5, 3, vax(counter));
    plot(xdata, data_wheel(:, trial));set(gca, 'XLim', [0 30]);
    if counter == 1
        title('velocity');
    end
end
if saveOn
    saveas(gcf, fullfile(savepath, [figName '.fig']));
    saveas(gcf, fullfile(savepath, [figName '.jpg']));    
end

%% try overlaying everything
xdata = TE.Wheel.xData;
mult = 0;
trials = [1:2];
trials = trials + 2 * mult;
data_wheel = TE.Wheel.data.V';
figName = ['pickTrials_combined_' num2str(trials(1)) ' +'];
ensureFigure('pickTrials', 1);

for counter = 1:2
    trial = trials(counter);
    subplot(2, 1, counter); 
    plot(xdata, data_chat(:,trial), 'g'); hold on;
    plot(xdata, data_dat(:,trial), 'r'); 
%     plot(xdata, nanzscore(data_pupil(:, trial)), 'b');
    plot(xdata, nanzscore(data_wheel(:, trial)), 'k');
    
    set(gca, 'XLim', [0 30], 'YLim', [-5 10]); ylabel(['trial ' num2str(trial)]); 
end
if saveOn
    saveas(gcf, fullfile(savepath, [figName '.fig']));
    saveas(gcf, fullfile(savepath, [figName '.jpg']));    
end
%% final example
trial = 34;
figName = ['coherence_ex_trial_' num2str(trial)];
ensureFigure(figName, 1);
plot(xdata, nanzscore(data_pupil(:, trial))-1, 'Color', [0.8 0.8 0.8], 'LineWidth', 2); hold on;
plot(xdata, data_chat(:, trial), 'g'); 
plot(xdata, data_dat(:, trial), 'r'); 
set(gca, 'YLim', [-3 7]);
legend({'ChAT', 'DAT', 'Pupil'}); legend('boxoff');
formatFigureGRC;
if saveOn
    saveas(gcf, fullfile(savepath, [figName '.fig']));
    saveas(gcf, fullfile(savepath, [figName '.jpg']));    
    saveas(gcf, fullfile(savepath, [figName '.epsc'])); 
end
