load(fullfile('Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk', 'DB.mat'));
% compile photometry averages, make grand average
fdField = 'ZS';
tau = 1.5; % time constant of exponential decay kernel
kDuration = 2.5; % duration of kernel
Fs = 20;
kt = (1:kDuration*20)*(1/Fs) - (1/Fs); 
k = exp(-1 * (1/tau) * kt);
% deconvolve data
epsilon = 0.1;

grand = struct();
figs = zeros(length(DB.animals), 2);
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadExperiment(DB, animal);
    if ~success
        disp('wtf');
        continue
    end    
    disp(animal);
    for channel = 1:2
        TE.Photometry.data(channel).([fdField 'deconv']) = bpDeconv(TE.Photometry.data(channel).(fdField), k, epsilon, 'none');
        avgData = phAverageFromTE(TE, {rewardTrials & csPlusTrials & hitTrials, uncuedReward}, channel,...
            'FluorDataField', [fdField 'deconv'], 'window', [-5, 2], 'zeroTimes', TE.Us); %high value, reward
            grand(channel).reward.components(:,:,counter) = avgData.Avg;
            grand(channel).reward.xData = avgData.xData(1,:);
%         
%         if channel == 2
%             subplot(3,3,counter)
%             plot(avgData.xData(1,:), avgData.Avg');
%         end

        avgData = phAverageFromTE(TE, {csPlusTrials & hitTrials, csPlusTrials & missTrials, csMinusTrials & CRTrials}, channel,...
            'FluorDataField', [fdField 'deconv'], 'window', [-2, 2], 'zeroTimes', TE.Cue); %high value, reward
            grand(channel).cue.components(:,:,counter) = avgData.Avg;
            grand(channel).cue.xData = avgData.xData(1,:);         
            
        % plot them
        figs(counter, channel) = ensureFigure(sprintf('deconvCompare_%s_ch%d', animal, channel), 1);
        subplot(1,2,1); 
        [~,~,ph] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, channel, 'FluorDataField', [fdField 'deconv'], 'window', [-4, 5]); set(gca, 'XLim', [-1 4]);
%         legend(ph, 'cuedReward', 'uncuedReward', 'Location','northwest', 'Box', 'off');
        xlabel('time from cue (s)'); ylabel('ZS'); title('raw'); textBox(animal);
        subplot(1,2,2); 
        [~,~,ph] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, channel, 'FluorDataField', [fdField 'deconv'], 'window', [-4, 5]); set(gca, 'XLim', [-1 4]);       
        legend(ph, {'cuedReward', 'uncuedReward'}, 'Location', 'northwest', 'Box', 'off')
        xlabel('time from cue (s)'); title('deconvolved');        
    end
end

% savepath = 'Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk\analysis\deconvolution_dev\';
% pdfname1 = fullfile(savepath,sprintf('deconv_tau=%d_ch1_eps0.pdf', tau));
% pdfname2 = fullfile(savepath,sprintf('deconv_tau=%d_ch2_eps0.pdf', tau));
% for counter = 1:size(figs, 1)
%     if counter == 1
%         export_fig(figs(counter, 1),pdfname1);  % write to pdf
%         export_fig(figs(counter, 2),pdfname2);  % write to pdf
%     else
%         export_fig(figs(counter, 1),'-append',pdfname1);  % write to pdf
%         export_fig(figs(counter, 2),'-append',pdfname2);  % write to pdf        
%     end
%         
% end

%%
ensureFigure('grandAverages_deconv', 1);
cueMap = [0 0 1; 0 1 1; 1 0 0];
rewMap = [0 0 1; 0 0 0];
subplot(2,2,1); title('Cue');
[~, hp] = boundedline(grand(1).cue.xData, squeeze(mean(grand(1).cue.components, 3))', permute(...
    std(grand(1).cue.components, 0,  3) / sqrt(size(grand(1).cue.components, 3)), [2 3 1]), 'cmap', cueMap);
legend(hp, {'CS+, hit', 'CS+, miss', 'CS-, CR'}, 'Location', 'northwest', 'Box', 'off');
ylabel('BF deconv');

subplot(2,2,3)
[~, hp] = boundedline(grand(2).cue.xData, squeeze(mean(grand(2).cue.components, 3))', permute(...
    std(grand(2).cue.components, 0,  3) / sqrt(size(grand(2).cue.components, 3)), [2 3 1]), 'cmap', cueMap);
legend(hp, {'CS+, hit', 'CS+, miss', 'CS-, CR'}, 'Location', 'northwest', 'Box', 'off');
xlabel('time from cue (s)'); ylabel('VTA deconv');

subplot(2,2,2); title('Reward');
[~, hp] = boundedline(grand(1).reward.xData, squeeze(mean(grand(1).reward.components, 3))', permute(...
    std(grand(1).reward.components, 0,  3) / sqrt(size(grand(1).reward.components, 3)), [2 3 1]), 'cmap', rewMap);
legend(hp, {'CS+, hit', 'uncued'}, 'Location', 'northwest', 'Box', 'off');

subplot(2,2,4)
[~, hp] = boundedline(grand(2).reward.xData, squeeze(mean(grand(2).reward.components, 3))', permute(...
    std(grand(2).reward.components, 0,  3) / sqrt(size(grand(2).reward.components, 3)), [2 3 1]), 'cmap', rewMap);
legend(hp, {'CS+, hit', 'uncued'}, 'Location', 'northwest', 'Box', 'off');
xlabel('time from reward (s)');