%% calculate latency, jitter, and reliability of uncued Us responses
% bl_window = [-1 0];
% response_window = [0 2]; % must start at 0
% PhotometryField = 'Photometry';
% Fs = TE.(PhotometryField).sampleRate;
% na = length(DB.animals);
% us_pooled = struct(...
%     'rew_latency', zeros(na, 2),...
%     'rew_jitter_std', zeros(na, 2),...
%     'rew_jitter_avg', zeros(na, 2),...
%     'rew_jitter_tt50', [],... % cell array to hold distributions
%     'rew_jitter_delta', [],... % cell array to hold distributions
%     'puff_latency', zeros(na, 2),...
%     'puff_jitter', zeros(na, 2),...
%     'puff_jitter_tt50', [],... % cell array to hold distributions    
%     'puff_jitter_delta', []... % cell array to hold distributions    
%     );
% for counter = 1:length(DB.animals)
%     animal = DB.animals{counter}
%     dbLoadAnimal(DB, animal);
%     theseTrials = find(rewardTrials & uncuedTrials & (TE.licks_us.rate > minRewardLickRate));    
%     [tt50, deltaMax] = deal(NaN(length(theseTrials), 2));
%     for channel = channels
%         [responseData, xData] = phAlignedWindow(TE, theseTrials, channel, 'zeroTimes', TE.Us, 'window', response_window, 'FluorDataField', 'ZS');
%         [blData, blx] = phAlignedWindow(TE, theseTrials, channel, 'zeroTimes', TE.Us, 'window', bl_window, 'FluorDataField', 'ZS');
%         [M, I] = max(responseData(:,zeroPoint:end), [], 2);
%         bl = mean(blData, 2);
%         deltaMax(:,channel) = M - bl;
%         level50 = bl + deltaMax(:,channel) * 0.5;
%         rew_avg = nanmean(responseData);        
%         [~, iAvg] = max(rew_avg);
%         us_pooled.rew_latency(counter, channel) = bpPnt2x(iAvg, Fs, 0);
%         % find time to 50% per trial and magnitude of trial response
% 
%         thresh = responseData >= level50;
%         crossings = [zeros(size(thresh,1), 1) thresh(:,2:end) - thresh(:,1:end-1)];
%         crossings = crossings > 0; % only want upward crossings
%         [row, col] = find(thresh);
%         pointMatrix = NaN(size(thresh));
%         pointMatrix(sub2ind(size(thresh), row,col)) = col;
%         valid = isfinite(min(pointMatrix,[],2)); % only include trials where crossing is found
%         theseCrossings = min(pointMatrix(valid,:),[],2);
%         tt50(valid,channel) = xData(theseCrossings);        
% 
% %         ensureFigure('test', 1);
% %         trials = randperm(length(theseTrials), 16);
% %         for counter = 1:16
% %             trial = trials(counter);
% %             subplot(4,4,counter); hold on;
% %         %     plot(xData, [responseData(trial,:)' zeros(numel(xData), 1) + level50(trial)]);
% %             plot(xData, responseData(trial,:));
% %             plot(blx, blData(trial,:));
% %             scatter(tt50(trial), level50(trial));
% %         end    
%     end       
%     us_pooled.rew_jitter_tt50{counter} = tt50;
%     us_pooled.rew_jitter_std(counter, :) = nanstd(tt50);
%     us_pooled.rew_jitter_avg(counter, :) = nanmean(tt50);
%     us_pooled.rew_jitter_delta{counter} = deltaMax;    
% end


%%
ensureFigure('test', 1);
trials = randperm(length(theseTrials), 16);
for counter = 1:16
    trial = trials(counter);
    subplot(4,4,counter); hold on;
%     plot(xData, [responseData(trial,:)' zeros(numel(xData), 1) + level50(trial)]);
    plot(xData, responseData(trial,:), '-');
    plot(blx, blData(trial,:));
    scatter(tt50(trial), level50(trial));
end
    