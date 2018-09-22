load(fullfile('Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk', 'DB.mat'));
% compile photometry averages, make grand average
fdField = 'ZS';
tau = 1; % time constant of exponential decay kernel
kDuration = 3; % duration of kernel
Fs = 20;
kt = (1:kDuration*20)*(1/Fs) - (1/Fs); % 2 second long kernel
k = exp(-1 * (1/tau) * kt);
% deconvolve data
epsilon = 0.1;


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
        TE.Photometry.data(channel).([fdField 'deconv']) = bpDeconv(TE.Photometry.data(channel).(fdField), k, epsilon);
        avgData = phAverageFromTE(TE, {rewardTrials & csPlusTrials & hitTrials, uncuedReward}, channel,...
            'FluorDataField', fdField, 'window', [1, 7]); %high value, reward
            grand(channel).reward.components(:,:,counter) = avgData.Avg;
            grand(channel).reward.xData = avgData.xData(1,:); 
%         
%         if channel == 2
%             subplot(3,3,counter)
%             plot(avgData.xData(1,:), avgData.Avg');
%         end

        avgData = phAverageFromTE(TE, {csPlusTrials & hitTrials, csPlusTrials & missTrials, csMinusTrials & CRTrials}, channel,...
            'FluorDataField', fdField, 'window', [-3, 7]); %high value, reward
            grand(channel).cue.components(:,:,counter) = avgData.Avg;
            grand(channel).cue.xData = avgData.xData(1,:);         
            
        % plot them
        figs(counter, channel) = ensureFigure(sprintf('deconvCompare_%s_ch%d', animal, channel), 1);
        subplot(1,2,1); 
        [~,~,ph] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, channel, 'FluorDataField', fdField, 'window', [-4, 5]); set(gca, 'XLim', [-1 4]);
%         legend(ph, 'cuedReward', 'uncuedReward', 'Location','northwest', 'Box', 'off');
        xlabel('time from cue (s)'); ylabel('ZS'); title('raw'); textBox(animal);
        subplot(1,2,2); 
        [~,~,ph] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, channel, 'FluorDataField', [fdField 'deconv'], 'window', [-4, 5]); set(gca, 'XLim', [-1 4]);       
        legend(ph, {'cuedReward', 'uncuedReward'}, 'Location', 'northwest', 'Box', 'off')
        xlabel('time from cue (s)'); title('deconvolved');
    end
end

savepath = 'Z:\SummaryAnalyses\LNL_Odor_v2_noPunish_whisk\analysis\deconvolution_dev\';
pdfname1 = fullfile(savepath,sprintf('deconv_tau=%d_ch1_eps0.pdf', tau));
pdfname2 = fullfile(savepath,sprintf('deconv_tau=%d_ch2_eps0.pdf', tau));
for counter = 1:size(figs, 1)
    if counter == 1
        export_fig(figs(counter, 1),pdfname1);  % write to pdf
        export_fig(figs(counter, 2),pdfname2);  % write to pdf
    else
        export_fig(figs(counter, 1),'-append',pdfname1);  % write to pdf
        export_fig(figs(counter, 2),'-append',pdfname2);  % write to pdf        
    end
        
end


% %% test deconvolution 
%     animal = 'DC_56';
%     dbLoadExperiment(DB, animal);
%     tau = 2; % time constant of exponential decay kernel
%     kDuration = 3; % duration of kernel
%     Fs = 20;
%     kt = (1:kDuration*20)*(1/Fs) - (1/Fs); % 2 second long kernel
%     k = exp(-1 * (1/tau) * kt);
%     % deconvolve data
%     epsilon = 0.1;
%     TE.Photometry.data(1).dFFdeconv = bpDeconv(TE.Photometry.data(1).dFF, k, epsilon);
%     TE.Photometry.data(2).dFFdeconv = bpDeconv(TE.Photometry.data(2).dFF, k, epsilon);    
%     ensureFigure('test_deconv', 1);
%     subplot(2,1,1);
%     phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, 1, 'FluorDataField', 'dFF', 'window', [-4, 7]);
%     subplot(2,1,2);
%     phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials, uncuedReward}, 1, 'FluorDataField', 'dFFdeconv', 'window', [-4, 7]);