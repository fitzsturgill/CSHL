


%% aligned to reward/beep
window = [-4 8];
channels = TE.Photometry.data.ch;

if ismember(channels, 2)    
    [rewards_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).raw, TE.Photometry.startTime, 20, TE.Reward, window);
end
[rewards_chat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(1).raw, TE.Photometry.startTime, 20, TE.Reward, window);

% local dFF
if ismember(channels, 2)
    bl_dat = nanmean(rewards_dat(:,1:40), 2);
    rewards_dat = bsxfun(@minus, rewards_dat, bl_dat);
    rewards_dat = bsxfun(@rdivide, rewards_dat, bl_dat);
end
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
if ismember(channels, 2)
    rewards_dat_sorted = rewards_dat(I, :);
end
rewards_chat_sorted = rewards_chat(I, :);

ensureFigure('random_rewards', 1); 
% subplot(3,2,1); imshow(rewards_chat_sorted, [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); colormap('jet'); title('chat'); ylabel('sorted');
% subplot(3,2,2); imshow(rewards_dat_sorted, [min(min(rewards_dat_sorted)), max(max(rewards_dat_sorted))]); colormap('jet'); title('dat');
subplot(3,2,1); image(rewards_chat, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); colormap('jet');  title('ChAT');
if ismember(channels, 2)
    subplot(3,2,2); image(rewards_dat, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', [min(min(rewards_dat_sorted)), max(max(rewards_dat_sorted))]); colormap('jet');    title('DAT');
end
xdata = linspace(window(1), window(2), size(rewards_chat, 2));
subplot(3,2,3); plot(xdata, nanmean(rewards_chat));
if ismember(channels, 2)
    subplot(3,2,4); plot(xdata, nanmean(rewards_dat));
end
subplot(3,2,5); triggeredEventRasterFromTE(TE, 'Port1In', TE.Reward);

if saveOn
    saveas(gcf, fullfile(savepath, 'random_rewards.fig'));
    saveas(gcf, fullfile(savepath, 'random_rewards.jpg'));
end


%% pupil data
%  [wheelY_new, wheelTimes_new] = resample(wheelY, wheelTimes, 20, 'linear');

TE = addPupilometryToTE(TE, 'duration', 30, 'zeroField', 'Baseline', 'startField', 'Baseline', 'frameRate', 60, 'frameRateNew', 20);


%% coherence
data_chat = TE.Photometry.data(1).raw';
% data_chat = data_chat(:,2:end) - data_chat(:,1:end-1); % whiten
% data_chat = nanzscore2(data_chat); % standardize
data_chat = nanzscore(data_chat); % standardize
if ismember(channels, 2)
    data_dat = TE.Photometry.data(2).raw';
    % data_dat = data_dat(:,2:end) - data_dat(:,1:end-1); % whiten
    % data_dat = nanzscore2(data_dat); % standardize
    data_dat = nanzscore(data_dat); % standardize
end

params.Fs = 20;
params.trialave = 1;
params.err = [2 0.05];
params.tapers = [3 5];
params.pad = 1;


if ismember(channels, 2)
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
end


%% coherence with pupil

data_pupil = TE.pupil.pupDiameter';
% data_pupil = data_pupil(:,2:end) - data_pupil(:,1:end-1); % whiten
% data_pupil = nanzscore2(data_pupil);
data_pupil = nanzscore(data_pupil);

validTrials = find(sum(isnan(data_pupil)) == 0);

params.Fs = 20;
params.trialave = 1;
params.err = [2 0.1];
params.tapers = [3 5];



ensureFigure('coherence_pupil', 1);

[C,phi,S12,S1,S2,f,confC, phistd, Cerr] = coherencyc(data_chat(:,validTrials), data_pupil(:,validTrials), params);
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


if ismember(channels, 2)    
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
end
if saveOn
    saveas(gcf, fullfile(savepath, 'coherence_pupil.fig'));
    saveas(gcf, fullfile(savepath, 'coherence_pupil.jpg'));
end
%%

ensureFigure('ch1_vs_pupil_examples', 1);
trial = 4;

subplot(2,1,1);
plot(TE.pupil.xData, TE.pupil.pupDiameter(trial,:)); ylabel('pupil');
subplot(2,1,2);
plot(TE.Photometry.xData, TE.Photometry.data.dFF(trial, :)); ylabel('GCaMP'); xlabel('time (s)');

%%

ensureFigure('ch1_vs_pupil_examples', 1);
trial = 4;

rawTimes = sessions(1).SessionData.RawEvents.Trial{trial}.Events.Port3In;
 wheelTimes = rawTimes(rawTimes < 30);
 wheelY = ones(size(wheelTimes));
 wheelTimes = [0 wheelTimes 30];
 wheelY = cumsum([0 wheelY 0]);
 
 [wheelY_new, wheelTimes_new] = resample(wheelY, wheelTimes, 20, 'linear');
 
 
 
D = 12.7; % diameter of wheel in cm
pr = 200; % pulses/rotation 
dpp = pi * 12.7 / 200;

wheelTimes = rawTimes(rawTimes < 30);
edges = 0:1/20:30;
% wheel_X = cumsum(histcounts(rawTimes, edges) * dpp);
wheel_X = smooth(cumsum(histcounts(rawTimes, edges)), 20);  
wheel_t = 0:1/20:30-1/20;
  
subplot(2,1,1);
plot(wheel_t, gradient(wheel_X)); ylabel('pupil');
subplot(2,1,2);
plot(TE.Photometry.xData, TE.Photometry.data.dFF(trial, :)); ylabel('GCaMP'); xlabel('time (s)');




