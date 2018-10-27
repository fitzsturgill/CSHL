% reversals_noPunish_exampleMouse

% using... DC_56


DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);
smoothWindow = 1;
saveOn = 1;
animal = 'DC_56';

photometryField = 'Photometry';
fdField = 'ZS';

success = dbLoadAnimal(DB, animal); % load TE and trial lookups

%% plot rasters for an example session/reversal

    
    saveName = sprintf('example_allBehavior_csPlus_%s', animal); 
    ensureFigure(saveName, 1);
    
    csPlusExampleTrials = csPlusTrials & rewardTrials & ismember(TE.sessionIndex, [2]);
    reversals = find(diff(TE.BlockNumber(csPlusExampleTrials, :))) + 1;
    sessionChanges = find(diff(TE.sessionIndex(csPlusExampleTrials, :))) + 1;

    
    climfactor = 2;
    subplot(1,5,1);       
    CData = TE.pupil.pupDiameterNorm(csPlusExampleTrials, :) - nanmean(TE.pupil.pupDiameterNorm(csPlusExampleTrials, 1:bpX2pnt(0,20,-4)), 2);
    image(CData, 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(CData(:)) - std(CData(:), 'omitnan') * climfactor, nanmean(CData(:)) + std(CData(:), 'omitnan') * climfactor]); 
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('pupil');  
    ylabel('trial number');
    
    
    subplot(1,5,2);
%     CData = TE.Whisk.whiskNorm(csPlusExampleTrials, :) - nanmean(TE.Whisk.whiskNorm(csPlusExampleTrials, 1:bpX2pnt(0,20,-4)), 2);
    imagesc(TE.Whisk.whiskNorm(csPlusExampleTrials, :), 'XData', [-4 7], [0 2])   
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'w', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('whisking');
    
    subplot(1,5,3);
    [~, lh] = eventRasterFromTE(TE, csPlusExampleTrials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording');
    line(repmat([-4; 7], 1, length(sessionChanges)), [sessionChanges'; sessionChanges'], 'Parent', gca, 'Color', 'k', 'LineWidth', 1); % reversal lines    
    line(repmat([-3; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines       
    set(gca, 'XLim', [-4 7]);
    set(lh, 'LineWidth', 0.3, 'Color', [0 0 0]);
    title('licking');
    xlabel('Time from odor (s)');
    climfactor = 3;  
    subplot(1,5,4); phRasterFromTE(TE, csPlusExampleTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('ChAT'); 
    
    subplot(1,5,5); phRasterFromTE(TE, csPlusExampleTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField); % 'CLimFactor', CLimFactor,
    line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 1); % reversal lines    
    title('DAT');    
    axs = findobj(gcf, 'Type', 'axes');    
    set(axs(1:end -1), 'YTick', []);
    
    formatFigurePublish('size', [4 1.5]);
    if saveOn 
        export_fig(fullfile(savepath, saveName), '-eps');
    end
    
%% averages median split on 

fdField = 'ZS';
saveName = sprintf('%s_phAvgs_%s', animal, fdField);  
h=ensureFigure(saveName, 1); 
mcLandscapeFigSetup(h);



subplot(2, 2, 1);
[ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, neutralTrials, uncuedReward}, 1,...
    'FluorDataField', fdField, 'window', [-4, 7], 'linespec', {'b','k','c'}); %high value, reward
legend(hl, {'rew', 'neutral', 'uncued rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
ylabel(sprintf('BF %s', fdField)); 



subplot(2, 2, 2);
[ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, neutralTrials, uncuedReward}, 2,...
    'FluorDataField', fdField, 'window', [-4, 7], 'linespec', {'b','k','c'}); %high value, reward
legend(hl, {'rew', 'neutral', 'uncued rew'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
ylabel(sprintf('VTA %s', fdField)); xlabel('time from cue (s)'); %set(gca, 'YLim', ylim);


subplot(2, 2, 3);
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 1,...
'FluorDataField', fdField, 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
title('CS+, outcomes'); set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);


subplot(2, 2, 4);
[ha, hl] = phPlotAverageFromTE(TE, {csPlusTrials & rewardTrials & hitTrials, csPlusTrials & rewardTrials & missTrials}, 2,...
    'FluorDataField', fdField, 'window', [-3, 7], 'linespec', {'c', 'm'}); %high value, reward
legend(hl, {'hit', 'miss'}, 'Location', 'southwest', 'FontSize', 12); legend('boxoff');
xlabel('time from cue (s)');     set(gca, 'XLim', [-3, 7]);%set(gca, 'YLim', ylim);

%% median split and plot histograms, not super convincing
bottom = (TE.licks_cs.rate == 0);
top = (TE.licks_cs.rate > 0);
bottom_ch1 = cum(TE.phPeakPercentile_us(1).data(bottom & rewardTrials));
top_ch1 = cum(TE.phPeakPercentile_us(1).data(top & rewardTrials));
bottom_ch2 = cum(TE.phPeakPercentile_us(2).data(bottom & rewardTrials));
top_ch2 = cum(TE.phPeakPercentile_us(2).data(top & rewardTrials));

ensureFigure('csLicks_medianSplit', 1);
subplot(1,2,1);
plot(bottom_ch1.sorted, bottom_ch1.index, ':', 'Color', mycolors('chat')); hold on;
plot(top_ch1.sorted, top_ch1.index, '-', 'Color', mycolors('chat')); 
plot(bottom_ch2.sorted, bottom_ch2.index, ':', 'Color', mycolors('dat')); 
plot(top_ch2.sorted, top_ch2.index, '-', 'Color', mycolors('dat')); 

subplot(1,2,2);
% scatter(TE.licks_cs.count(rewardTrials), TE.phPeakMean_us(2).data(rewardTrials), '.');
lickBins = [0:5 10] - 0.5;
[ch1_means, ch1_sem, ~] = binnedMeansXY(TE.licks_cs.count(rewardTrials), TE.phPeakMean_us(1).data(rewardTrials), lickBins);
[ch2_means, ch2_sem, binCenters] = binnedMeansXY(TE.licks_cs.count(rewardTrials), TE.phPeakMean_us(2).data(rewardTrials), lickBins);
errorbar(binCenters, ch1_means, ch1_sem, 'Color', mycolors('chat')); hold on;
errorbar(binCenters, ch2_means, ch2_sem, 'Color', mycolors('dat'));
% if saveOn
%     saveas(gcf, fullfile(savepath, [saveName '.fig']));
%     saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
% end    




%% Reversal Averages


% load the reversals
load(fullfile(DB.path, 'pooled', sprintf('RE_%s.mat', animal)));    
%% find good reversals
% best reversal is reversal #1
smoothWindow = 3;
% plot individual reversals
nRev = size(RE.csPlus.phPeakMean_cs_ch1.before, 1);
ensureFigure('findGoodRevs', 1);
xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
% bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);
% ceilIx(1) = nearest(RE.csMinus.trialsBefore, -30); ceilIx(2) = nearest(RE.csPlus.trialsBefore, 0);
for counter = 1:nRev
    subplot(3,4,counter); hold on;
    plot(xData, smoothdata([RE.csMinus.licks_cs.before(counter,:) RE.csPlus.licks_cs.after(counter, :)], 2, 'movmean', smoothWindow),'k'); hold on;
%     plot(xData, smoothdata([RE.csMinus.phPeakMean_cs_ch1.before(counter,:) RE.csPlus.phPeakMean_cs_ch1.after(counter, :)], 2, 'movmean', smoothWindow), 'g');
%     plot(xData, smoothdata([RE.csMinus.phPeakMean_cs_ch2.before(counter,:) RE.csPlus.phPeakMean_cs_ch2.after(counter, :)], 2, 'movmean', smoothWindow), 'r');
%     plot(xData, smoothdata([RE.csMinus.licks_cs.before(counter,:) RE.csPlus.licks_cs.after(counter, :)], 2, 'movmean', 1),'k.'); hold on;
    plot(xData, smoothdata([RE.csMinus.phPeakMean_cs_ch1.before(counter,:) RE.csPlus.phPeakMean_cs_ch1.after(counter, :)], 2, 'movmean', 1), 'g.');
    plot(xData, smoothdata([RE.csMinus.phPeakMean_cs_ch2.before(counter,:) RE.csPlus.phPeakMean_cs_ch2.after(counter, :)], 2, 'movmean', 1), 'r.');
    textBox(int2str(counter));
end
%% plot reversal # 1 and fit a weibull function to it
smoothWindow = 3;
% plot individual reversals
nRev = size(RE.csPlus.phPeakMean_cs_ch1.before, 1);
rev = 1;
saveName = 'exampleReversal_newCsPlus';
ensureFigure(saveName, 1);
ch1Data = [RE.csMinus.phPeakMean_cs_ch1.before(rev,:) RE.csPlus.phPeakMean_cs_ch1.after(rev, :)];
ch2Data = [RE.csMinus.phPeakMean_cs_ch2.before(rev,:) RE.csPlus.phPeakMean_cs_ch2.after(rev, :)];
lickData = [RE.csMinus.licks_cs.before(rev,:) RE.csPlus.licks_cs.after(rev, :)];
validtrials = ~isnan(ch1Data);
ch1Data = ch1Data(validtrials);
ch2Data = ch2Data(validtrials);
lickData = lickData(validtrials);

xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter] + 1; % doesn't make sense that last trial of last reversal is trial 0 does it?
xData = xData(validtrials);
model = 'a * (1 - exp(-1 * (x/b)^c)) + d'; % weibull function, CDF form
fo = fitoptions('Method', 'NonlinearLeastSquares',...
    'Upper', [Inf  Inf Inf Inf],...
    'Lower', [0 0 0 -Inf],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
    'StartPoint', [1 round(length(ch1Data)/2) 1 0]...
    );
ft = fittype(model, 'options', fo);
%     xData = (0:length(toFit) - 1)';
[fitobject_ch1, ~, ~] = fit((0:(length(ch1Data) - 1))', ch1Data', ft, fo);
[fitobject_ch2, ~, ~] = fit((0:(length(ch1Data) - 1))', ch2Data', ft, fo);

axes; hold on;
hl = [];
hl(end + 1) = plot(xData, ch1Data, '.', 'Color', mycolors('chat'), 'MarkerSize', 4);
hl(end + 1) = plot(xData, ch2Data, '.', 'Color', mycolors('dat'), 'MarkerSize', 4);
plot(xData, fitobject_ch1.a .* (1 - exp(-1 * ((xData - xData(1)) ./ fitobject_ch1.b) .^fitobject_ch1.c)) + fitobject_ch1.d, ':', 'Color', mycolors('chat'));
plot(xData, fitobject_ch2.a .* (1 - exp(-1 * ((xData - xData(1)) ./ fitobject_ch2.b) .^fitobject_ch2.c)) + fitobject_ch2.d, ':', 'Color', mycolors('dat'));
set(gca, 'YLim', [-4 7]);    
yyaxis right; set(gca, 'YColor', [0 0 0]); ylabel('licks/s');
hl(end + 1) = plot(xData, smoothdata(lickData, 2, 'movmedian', 5),'k');
set(gca, 'YLim', [0 6], 'XLim', xData([1 end]));
yyaxis left;
xlabel('New Cs+ trials from reversal')
ylabel('\fontsize{8}Fluor. (\fontsize{12}\sigma\fontsize{8}-baseline)');
% legend({'ACh.', 'Dop.', 'Licks'}, 'Box', 'off', 'Location', 'best', 'Orientation', 'vertical');
formatFigurePublish('size', [2 1.1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end





%% plot the averages

channels = [1 2];
peakFieldCh1 = 'phPeakMean_cs_ch1';
peakFieldCh2 = 'phPeakMean_cs_ch2';    

saveName = [animal '_' strtok(strtok(peakFieldCh1, '_'), '_') '_phCue_revAvg'];

h=ensureFigure(saveName, 1);
mcLandscapeFigSetup(h);


% data fields from RE to plot
dataFields = {
    'licks_cs', 'k';...        
    peakFieldCh1, 'g';...
    peakFieldCh2, 'r';
    };


% CS MINUS -> CS PLUS
    % x data and indices for baseline (cs- baseline) and ceiling (cs+
    % baseline)
xData = [RE.csMinus.trialsBefore RE.csPlus.trialsAfter];
bl(1) = nearest(xData, -30);bl(2) = nearest(xData, 0);
ceilIx(1) = nearest(RE.csMinus.trialsBefore, -30); ceilIx(2) = nearest(RE.csPlus.trialsBefore, 0);

% scale data and plot
nFields = size(dataFields, 1);
for counter = 1:nFields
    dataField = dataFields{counter, 1};
    color = dataFields{counter, 2};
    revNorm = [RE.csMinus.(dataField).before RE.csPlus.(dataField).after];
    revNorm = smoothdata(revNorm, 2, 'movmean', smoothWindow, 'omitnan');
%         baseline = nanmean(RE.csMinus.(dataField).before(:,bl(1):bl(2)), 2); % subtract off Cs- prior to reversal
%     ceiling = percentile(RE.csPlus.(dataField).before(:,ceilIx(1):ceilIx(2)), 0.8, 2); % divide by Cs+ prior to reversal
%         revNorm = revNorm - baseline; % RIGHT NOW JUST NORMALIZE, NO BASELINE SUBTRACTION    
%     revNorm = revNorm ./ ceiling;
    subplot(2, 4, 1); hold on;
    boundedline(xData, nanmean(revNorm), nanSEM(revNorm), color);
    set(gca, 'XLim', [-50 50]);
    subplot(2,4,counter + 1);
    imagesc('XData', xData, 'CData', revNorm);
    set(gca, 'XLim', [-50 50]);
end
%     set(gca, 'XLim', [-50 50]);
%     xlabel('Trials of new CS+ odor from reversal'); 
%     ylabel('Cue response normalized'); title([strtok(peakFieldCh1, '_') ', avg ' num2str(nReversals) ' reversals'], 'Interpreter', 'none');

% CS PLUS -> CS MINUS
subplot(2,1,2);
    % x data and indices for baseline (cs- baseline) and ceiling (cs+
    % baseline)        
xData = [RE.csPlus.trialsBefore RE.csMinus.trialsAfter];
bl(1) = nearest(RE.csMinus.trialsBefore, -30); bl(2) = nearest(RE.csMinus.trialsBefore, 0);
ceilIx(1) = nearest(RE.csPlus.trialsBefore, -30); ceilIx(2) = nearest(RE.csPlus.trialsBefore, 0);

% scale data and plot
for counter = 1:size(dataFields, 1)
    dataField = dataFields{counter, 1};
    color = dataFields{counter, 2};
    revNorm = [RE.csPlus.(dataField).before RE.csMinus.(dataField).after];
    revNorm = smoothdata(revNorm, 2, 'movmean', smoothWindow, 'omitnan');
%         baseline = nanmean(RE.csMinus.(dataField).before(:,bl(1):bl(2)), 2); % subtract off Cs- prior to reversal
    ceiling = percentile(RE.csPlus.(dataField).before(:,ceilIx(1):ceilIx(2)), 0.8, 2); % divide by Cs+ prior to reversal
%         revNorm = revNorm - baseline; % RIGHT NOW JUST NORMALIZE, NO BASELINE SUBTRACTION         
    revNorm = revNorm ./ ceiling;
    subplot(2,4,5); hold on;
    boundedline(xData, nanmean(revNorm), nanSEM(revNorm), color);
    set(gca, 'XLim', [-50 50]);
    subplot(2,4,4 + counter + 1);
    imagesc('XData', xData, 'CData', revNorm);
    set(gca, 'XLim', [-50 50]);
end
%     set(gca, 'XLim', [-50 50]);
%     xlabel('Trials of new CS- odor from reversal'); 
%     ylabel('Cue response normalized'); title([strtok(peakFieldCh1, '_') ', avg ' num2str(nReversals) ' reversals'], 'Interpreter', 'none');
    
% 
%     if saveOn
%         saveas(gcf, fullfile(savepath, [saveName '.fig']));
%         saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
%     end