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
rawTimes = sessions(1).SessionData.RawEvents.Trial{5}.Events.Port3In;
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

[rewards_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).raw, TE.Photometry.startTime, 20, TE.Reward, [-3 3]);
rewards_chat = extractDataByTimeStamps(TE.Photometry.data(1).raw, TE.Photometry.startTime, 20, TE.Reward, [-3 3]);

% local dFF
bl_dat = nanmean(rewards_dat(:,1:60), 2);
rewards_dat = bsxfun(@minus, rewards_dat, bl_dat);
rewards_dat = bsxfun(@rdivide, rewards_dat, bl_dat);
bl_chat = nanmean(rewards_chat(:,1:60), 2);
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

ensureFigure('rewards_sorted', 1); 
subplot(3,2,1); imshow(rewards_chat_sorted, [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); colormap('jet'); title('chat'); ylabel('sorted');
subplot(3,2,2); imshow(rewards_dat_sorted, [min(min(rewards_dat_sorted)), max(max(rewards_dat_sorted))]); colormap('jet'); title('dat');
subplot(3,2,3); imshow(rewards_chat, [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); colormap('jet'); ylabel('unsorted');
subplot(3,2,4); imshow(rewards_dat, [min(min(rewards_dat_sorted)), max(max(rewards_dat_sorted))]); colormap('jet');
subplot(3,2,5); plot(nanmean(rewards_chat));
subplot(3,2,6); plot(nanmean(rewards_dat));



%%

[rewards_wheel, ts, tn] = extractDataByTimeStamps(TE.Wheel.data.V, TE.Wheel.startTime, 20, TE.Reward, [-1 1]);
ensureFigure('rewards_wheel', 1);
plot(nanmean(rewards_wheel));

%% plot individual trials with reward times annotated

channel = 2;
ensureFigure(['annotated_' num2str(channel)], 1);
for trial = 1:18

    subplot(6,3,trial);
    trial = trial + 18 * 0;
    ydata = TE.Photometry.data(channel).raw(trial, :);    
    plot(TE.Photometry.xData, ydata, 'k'); hold on;
    tsx = repmat(TE.Reward{trial}(:,1), 1, 2)';
    tsy = [repmat(min(min(ydata)), 1, size(tsx, 2)); repmat(max(max(ydata)), 1, size(tsx, 2))];    
    plot(tsx, tsy, 'r');
end











 