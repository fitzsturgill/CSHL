window = [-2 2];
% [rewards_dat, ts, tn] = extractDataByTimeStamps(TE.Photometry.data(2).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2]);
rewards_chat = extractDataByTimeStamps(TE.Photometry.data(1).raw, TE.Photometry.startTime, 20, TE.Reward, [-2 2]);

% local dFF
% bl_dat = nanmean(rewards_dat(:,1:40), 2);
% rewards_dat = bsxfun(@minus, rewards_dat, bl_dat);
% rewards_dat = bsxfun(@rdivide, rewards_dat, bl_dat);
bl_chat = nanmean(rewards_chat(:,1:40), 2);
rewards_chat = bsxfun(@minus, rewards_chat, bl_chat);
rewards_chat = bsxfun(@rdivide, rewards_chat, bl_chat);
% 

ts_abs = TE.trialStartTimeStamp;

iri_pre = [Inf; diff(ts_abs)];
iri_post = [diff(ts_abs); Inf];

[~, I] = sort(iri_pre);
iri_pre_sorted = iri_pre(I);
% rewards_dat_sorted = rewards_dat(I, :);
rewards_chat_sorted = rewards_chat(I, :);

ensureFigure('random_rewards', 1); 
% subplot(3,2,1); imshow(rewards_chat_sorted, [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); colormap('jet'); title('chat'); ylabel('sorted');
% subplot(3,2,2); imshow(rewards_dat_sorted, [min(min(rewards_dat_sorted)), max(max(rewards_dat_sorted))]); colormap('jet'); title('dat');
subplot(3,2,1); image(rewards_chat, 'XData', window, 'CDataMapping', 'Scaled'); set(gca, 'CLim', [min(min(rewards_chat_sorted)), max(max(rewards_chat_sorted))]); colormap('jet');  title('ChAT');
% set(gca, 'CLim', [-0.004 0.004]);
% subplot(3,2,2); image(rewards_dat, 'XData', window,  'CDataMapping', 'Scaled'); set(gca, 'CLim', [min(min(rewards_dat_sorted)), max(max(rewards_dat_sorted))]); colormap('jet');    title('DAT');
xdata = linspace(window(1), window(2), size(rewards_chat, 2));
subplot(3,2,3); plot(xdata, nanmean(rewards_chat));
% subplot(3,2,4); plot(xdata, nanmean(rewards_dat));
subplot(3,2,5); triggeredEventRasterFromTE(TE, 'Port1In', TE.Reward);

if saveOn
    saveas(gcf, fullfile(savepath, 'random_rewards.fig'));
    saveas(gcf, fullfile(savepath, 'random_rewards.jpg'));
end
