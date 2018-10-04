load(fullfile('Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk', 'DB.mat'));
% compile photometry averages, make grand average
fdField = 'ZS';
tau = 2; % time constant of exponential decay kernel
kDuration = 3; % duration of kernel
Fs = 20;
kt = (1:kDuration*20)*(1/Fs) - (1/Fs); 
k = exp(-1 * (1/tau) * kt);
k = k / trapz(k);
% deconvolve data
epsilon = 0.1;

grand = struct();
figs = [];
for counter = 1:length(DB.animals)
    animal = DB.animals{counter};
    success = dbLoadAnimal(DB, animal);
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
        figs(end + 1, channel) = ensureFigure(sprintf('deconvCompare_%s_ch%d', animal, channel), 1);
        subplot(1,2,1); 
        [~,~,ph] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, channel, 'FluorDataField', [fdField], 'window', [-4, 5]); set(gca, 'XLim', [-1 4]);
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

%% for lab meeting
animal = 'DC_46';
success = dbLoadAnimal(DB, animal);
figs = [];
fdField = 'dF';
TE.Photometry.data(1).([fdField 'deconv']) = bpDeconv(TE.Photometry.data(1).(fdField), k, 0, 'none');
TE.Photometry.data(2).([fdField 'deconv']) = bpDeconv(TE.Photometry.data(2).(fdField), k, 0, 'none');    
TE.Photometry.data(1).([fdField 'deconvReg']) = bpDeconv(TE.Photometry.data(1).(fdField), k, 0.2, 'none');
TE.Photometry.data(2).([fdField 'deconvReg']) = bpDeconv(TE.Photometry.data(2).(fdField), k, 0.2, 'none');
window = [-2 6];

% plot them
figs(end + 1) = ensureFigure('regularization_effect_ch1', 1);
subplot(2,2,1); 
[~,~,ph] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, 1, 'FluorDataField', [fdField], 'window', window, 'linespec', {'b', 'k'});
%         legend(ph, 'cuedReward', 'uncuedReward', 'Location','northwest', 'Box', 'off');
ylabel(fdField); title('raw'); textBox(animal);
subplot(2,2,2); 
[~,~,ph] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, 1, 'FluorDataField', [fdField 'deconv'], 'window', window, 'linespec', {'b', 'k'});
legend(ph, {'cuedReward', 'uncuedReward'}, 'Location', 'northwest', 'Box', 'off')
xlabel('time from cue (s)'); title('deconvolved');
subplot(2,2,3); 
[~,~,ph] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, 1, 'FluorDataField', [fdField 'deconvReg'], 'window', window, 'linespec', {'b', 'k'});   
xlabel('time from cue (s)'); title('deconvolved'); ylabel(fdField);
formatFigureTalk([4 2] * 2.5);

% kernel, raw and regularized
figs(end + 1) = ensureFigure('kernel', 1);
Lx2 = pow2(nextpow2(length(k)));
K = fft(k, Lx2); % transform kernel, kernel is constant so this only need be done once
Kreg = K + epsilon * mean(K.*conj(K))./conj(K);
Ksumm = epsilon * mean(K.*conj(K))./conj(K);
kreg = real(ifft(Kreg));
kf = 20 * (0:(Lx2/2))/Lx2;
subplot(1,2,1); hold on;
plot(kt, k, 'k', 'LineWidth', 2);
plot(kt, kreg(1:length(k)), 'g--', 'LineWidth', 2);
legend('raw', 'regularized', 'Box', 'off', 'Location', 'northeast');
title('Kernel, Time Domain'); xlabel('time (s)');
subplot(1,2,2); hold on;
plot(kf, K(1:length(kf)), 'k', 'LineWidth', 2);
plot(kf, Kreg(1:length(kf)), 'g', 'LineWidth', 2);
plot(kf, Ksumm(1:length(kf)), 'r', 'LineWidth', 2);
legend('raw', 'regularized', 'regularizer', 'Box', 'off', 'Location', 'northeast');
title('Kernel, Frequency Domain'); xlabel('Frequency (1/s)');
formatFigureTalk([4 2] * 2.5);


figs(end + 1) = ensureFigure('grandAverages_deconv', 1);
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
formatFigureTalk([4 3] * 2.5);