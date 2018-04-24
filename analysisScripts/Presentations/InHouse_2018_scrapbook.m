
inhouseSavepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\Presentations\2018_InHouse_Seminar';

%% sorted rasters snippet, in progress...
% ChAT_39 -> 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete_newPupil\ChAT_39'
    TE.firstLick = calcEventLatency(TE, 'Port1In', TE.Cue, TE.Us); % 3 seconds from start of odor to outcome, if there are no anticipatory licks, then call it 3 second latency to first lick
    CLimFactor = 1.5;
    figname = 'phRastersFromTE_reward_sorted';
    h=ensureFigure(figname, 1);
    xlim = [-5 4];
    trials = highValueTrials;
    subplot(1,2,1);    
    eventRasterFromTE(TE, trials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.firstLick);    
    addStimulusBar(gca, [-3 -2 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [-0.1 0.1 0], '', [0.3 0.3 0.3], 5);
    set(gca, 'XLim', xlim); xlabel('Time from reinforcement (s)'); set(gca, 'YTick', []);
    
    subplot(1,2,2);
    phRasterFromTE(TE, trials, 1, 'CLimFactor', CLimFactor, 'sortValues', TE.firstLick); hold on;
    sortValues = sort(TE.firstLick(trials));
    plot(sortValues + -3, 1:sum(trials), 'r', 'LineWidth', 2);
    set(gca, 'YTick', [], 'XLim', xlim);%, 'YTickLabel', {''});
    addStimulusBar(gca, [-3 -2 2], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [-0.1 0.1 2], '', [0.3 0.3 0.3], 5);
    xlabel('Time from reinforcement (s)');
%     title([TE.filename{1}(1:7) ': hival, reward'], 'Interpreter', 'none'); 
    
    formatFigureTalk([7 4]);

%     trials = lowValueTrials;
%     subplot(1,4,3);
%     eventRasterFromTE(TE, trials, 'Port1In', 'trialNumbering', 'consecutive',...
%         'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.firstLick);
%     
%     subplot(1,4,4);
%     phRasterFromTE(TE, trials, 1, 'CLimFactor', CLimFactor, 'sortValues', TE.firstLick); hold on;
%     sortValues = sort(TE.firstLick(trials));
%     plot(sortValues + -3, 1:sum(trials), 'r', 'LineWidth', 2);
% %     title([TE.filename{1}(1:7) ': hival, reward'], 'Interpreter', 'none'); 

    if saveOn
        saveas(gcf, fullfile(inhouseSavepath, figname), 'fig');
        saveas(gcf, fullfile(inhouseSavepath, figname), 'jpeg');   
        saveas(gcf, fullfile(inhouseSavepath, figname), 'emf');
        saveas(gcf, fullfile(inhouseSavepath, figname), 'epsc');           
    end
%%
load('Z:\SummaryAnalyses\CuedOutcome_Odor_Complete_newPupil\ChAT_39\TE.mat');
cuedOutcome_Conditions;

%%
% ChAT_39 -> 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete_newPupil\ChAT_39'
    saveName = 'all_behavior_highValue_withEye';
    ensureFigure(saveName, 1);
    subplot(1,4,1);
    
    imagesc(TE.pupil.eyeAreaNorm(highValueTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
%     set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
%     title('Velocity');
    subplot(1,4,2);
    image(TE.pupil.pupDiameterNorm(highValueTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('Pupil Diameter');    
    subplot(1,4,3); phRasterFromTE(TE, highValueTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
    ax = findobj(gcf, 'Type', 'Axes');
%     set(ax, 'YLim', [40 80]);
if saveOn
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));    
    disp('figure saved');
end

saveName = 'all_behavior_CsMinus';
    ensureFigure(saveName, 1);
%     reversals = find(diff(TE.BlockNumber(lowValueTrials, :))) + 1;
    subplot(1,4,1);
    title('area');
    imagesc(TE.pupil.eyeAreaNorm(lowValueTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
%     set(gca, 'CLim', [min(TE.Wheel.data.V(:)), max(TE.Wheel.data.V(:))]); 
%     set(gca, 'CLim', [mean(TE.Wheel.data.V(:)) - std(TE.Wheel.data.V(:)) * 2, mean(TE.Wheel.data.V(:)) + std(TE.Wheel.data.V(:)) * 2]); 
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    subplot(1,4,2);
    image(TE.pupil.pupDiameterNorm(lowValueTrials, :), 'XData', [-4 7], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [nanmean(TE.pupil.pupDiameterNorm(:)) - std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2, nanmean(TE.pupil.pupDiameterNorm(:)) + std(TE.pupil.pupDiameterNorm(:), 'omitnan') * 2]); 
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    colormap('parula');  
    title('Pupil Diameter');    
    subplot(1,4,3); phRasterFromTE(TE, lowValueTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
%     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
    title('ChAT'); xlabel('Time frome odor (s)');
%     subplot(1,4,4); phRasterFromTE(TE, lowValueTrials, 2, 'trialNumbering', 'consecutive', 'CLimFactor', 2); % 'CLimFactor', CLimFactor,
% %     line(repmat([-4; 7], 1, length(reversals)), [reversals'; reversals'], 'Parent', gca, 'Color', 'r', 'LineWidth', 2); % reversal lines    
%     title('DAT');    
if saveOn
    saveas(gcf, fullfile(inhouseSavepath, [saveName '.fig']));
    saveas(gcf, fullfile(inhouseSavepath, [saveName '.jpg']));    
    disp('figure saved');
end

%% cuedOutcome sort by ChAT
% ChAT_39 -> 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete_newPupil\ChAT_39'
    CLimFactor = 1.5;
    figname = 'cuedOutcome_rasters_ChATSorted';
    h=ensureFigure(figname, 1);
    xlim = [-4 4];
    trials = highValueTrials;
    trialsSorted = find(trials);
    sortValues = TE.phPeak_cs.data;
    [sorted index] = sort(TE.phPeak_cs.data(trials));
    trialsSorted = trialsSorted(index);
    
    subplot(1,3,3);
    image(TE.pupil.pupDiameterNorm(trialsSorted, :), 'XData', [-7 4], 'CDataMapping', 'Scaled');
    set(gca, 'CLim', [.5 1.5]); 
        set(gca, 'YTick', [], 'XLim', xlim);%, 'YTickLabel', {''});
    addStimulusBar(gca, [-3 -2 2], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [-0.1 0.1 2], '', [0.3 0.3 0.3], 5);
    
    subplot(1,3,2);    
    eventRasterFromTE(TE, trials, 'Port1In', 'trialNumbering', 'consecutive',...
        'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', sortValues);    
    addStimulusBar(gca, [-3 -2 0], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [-0.1 0.1 0], '', [0.3 0.3 0.3], 5);
    set(gca, 'XLim', xlim); xlabel('Time from reinforcement (s)'); set(gca, 'YTick', []);
    
    subplot(1,3,1);
    phRasterFromTE(TE, trials, 1, 'CLimFactor', CLimFactor, 'sortValues', sortValues); hold on;    
    set(gca, 'YTick', [], 'XLim', xlim);%, 'YTickLabel', {''});
    addStimulusBar(gca, [-3 -2 2], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [-0.1 0.1 2], '', [0.3 0.3 0.3], 5);
%     xlabel('Time from reinforcement (s)');
%     title([TE.filename{1}(1:7) ': hival, reward'], 'Interpreter', 'none'); 
    
    formatFigureTalk([7 4]);

%     trials = lowValueTrials;
%     subplot(1,4,3);
%     eventRasterFromTE(TE, trials, 'Port1In', 'trialNumbering', 'consecutive',...
%         'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.firstLick);
%     
%     subplot(1,4,4);
%     phRasterFromTE(TE, trials, 1, 'CLimFactor', CLimFactor, 'sortValues', TE.firstLick); hold on;
%     sortValues = sort(TE.firstLick(trials));
%     plot(sortValues + -3, 1:sum(trials), 'r', 'LineWidth', 2);
% %     title([TE.filename{1}(1:7) ': hival, reward'], 'Interpreter', 'none'); 

    if saveOn
        saveas(gcf, fullfile(inhouseSavepath, figname), 'fig');
        saveas(gcf, fullfile(inhouseSavepath, figname), 'jpeg');   
        saveas(gcf, fullfile(inhouseSavepath, figname), 'emf');
        saveas(gcf, fullfile(inhouseSavepath, figname), 'epsc');           
    end
    
    %% scatter plots
    % ChAT_39 -> 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete_newPupil\ChAT_39'
    figname = 'ChAT_predictive_scatter';
    ensureFigure(figname, 1);
    subplot(2,1,1);
    yData = TE.pupil_cs(highValueTrials);
    xData = TE.phPeak_cs_phasic.data(highValueTrials);
    include = ~isnan(yData);
    scatter(xData(include),yData(include)); hold on;
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData(include),yData(include), 'poly1', fo); 
    fph=plot(fob,'predfunc'); legend off;
    set(fph, 'LineWidth', 2);
    textBox(['R = ' num2str(corr(xData(include), yData(include)), 3)],[], [0.2 0.95], 16);
    set(gca, 'XLim', [-1 4], 'XTick', []);
%     xlabel('Cholinergic phasic (\fontsize{20}\sigma\fontsize{16}-baseline)');
    xlabel('');
    ylabel('Pupil diameter (norm.)');
    
    subplot(2,1,2);
    yData = TE.csLicks.rate(highValueTrials);
    xData = TE.phPeak_cs_phasic.data(highValueTrials);
    include = ~isnan(yData);
    scatter(xData(include),yData(include)); hold on;
    fo = fitoptions('poly1');%, 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);
    fob = fit(xData(include),yData(include), 'poly1', fo); 
    fph = plot(fob,'predfunc'); legend off;    
    set(fph, 'LineWidth', 2);
    set(gca, 'XLim', [-1 4]);
    xlabel('Cholinergic phasic (\fontsize{20}\sigma\fontsize{16}-baseline)');
    ylabel('Antic. Licks (1/s)');
    textBox(['R = ' num2str(corr(xData(include), yData(include)), 3)],[], [0.2 0.95], 16);
    formatFigureTalk([5 6]);
    
    if saveOn
        saveas(gcf, fullfile(inhouseSavepath, figname), 'fig');
        saveas(gcf, fullfile(inhouseSavepath, figname), 'jpeg');   
        saveas(gcf, fullfile(inhouseSavepath, figname), 'emf');
        saveas(gcf, fullfile(inhouseSavepath, figname), 'epsc');           
    end
    
    %% deconvolution to reveal outcome modulation (perhaps)
    uncuedPunishment = filterTE(TE, 'trialType', 8);
    kernel = phAverageFromTE(TE, uncuedPunishment, 1, 'window', [0 2], 'FluorDataField', 'ZS');
    data =  TE.Photometry.data.ZS;
    h = kernel.Avg; 
    TE.Photometry.data.deconv = deconv_Fourier(data, h, 0.2);
    
    
    
    figname = 'deconvolved';
    ensureFigure(figname, 1);
    subplot(1,2,1, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1, 'FluorDataField', 'ZS', 'window', [-5 2]); % reward, varying degrees of expectation
%     legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'northwest', 'FontSize', 16); legend('boxoff');
    title(''); ylabel('Z Score'); xlabel('time from reinforcement (s)');   
    set(gca, 'XLim', [-4 2], 'YLim', [-2 8]);
    addStimulusBar(gca, [-3 -2 8], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [-0.1 0.1 8], '', [0.3 0.3 0.3], 5);
    
    axes('Position', [0.3 0.7 0.15 0.2]); % 0.1300    0.1100    0.3347    0.8150
    phPlotAverageFromTE(TE, uncuedPunishment, 1, 'FluorDataField', 'ZS', 'window', [0 2], 'linespec', 'k');
    textBox('Kernel: punish', [], [0.6 1], 16);
    
    subplot(1,2,2, 'FontSize', 12, 'LineWidth', 1); [ha, hl] = phPlotAverageFromTE(TE, trialsByType([1 4 7]), 1, 'FluorDataField', 'deconv', 'window', [-5 1.9]); % reward, varying degrees of expectation
    addStimulusBar(gca, [-3 -2 8], '', [0.3 0.3 0.3], 5); addStimulusBar(gca, [-0.1 0.1 8], '', [0.3 0.3 0.3], 5);
    legend(hl, {'hival, rew', 'loval, rew', 'rew'}, 'Location', 'northwest', 'FontSize', 16); legend('boxoff');
    ylabel('Z Score'); xlabel('time from reinforcement (s)');   
    set(gca, 'XLim', [-4 2], 'YLim', [-2 8]);
    formatFigureTalk([10 5.5]);
    
    if saveOn
        saveas(gcf, fullfile(inhouseSavepath, figname), 'fig');
        saveas(gcf, fullfile(inhouseSavepath, figname), 'jpeg');   
        saveas(gcf, fullfile(inhouseSavepath, figname), 'emf');
        saveas(gcf, fullfile(inhouseSavepath, figname), 'epsc');           
    end
    
