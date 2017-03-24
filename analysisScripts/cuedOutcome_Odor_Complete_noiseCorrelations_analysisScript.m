% cuedOutcome_odor_complete --> noise correlation analysis

% assumes that you have already loaded the TE that exists in
% SummaryAnalyses
    %%
% savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\SummaryAnalyses\CuedOutcome_Odor_Complete';
% savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
% basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';
% basepath = uigetdir;
basepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete';
saveOn = 1;
%%
subjectName = TE.filename{1}(1:7);
disp(subjectName);
savepath = fullfile(basepath, subjectName);
ensureDirectory(savepath);
%%
cueLickWindow = [1 3]; % referenced from cue onset
usLickWindow = [0 2]; % referenced from us onset
cuePhWindow = [1 3];
usPhWindow = [0 0.5];

%% extract peak trial dFF responses to cues and reinforcement and lick counts
TE.phPeak_cs = bpCalcPeak_dFF(TE.Photometry, 1, cuePhWindow, TE.Cue, 'method', 'mean');
TE.phPeak_us = bpCalcPeak_dFF(TE.Photometry, 1, usPhWindow, TE.Us, 'method', 'mean');

TE.csLicks = countEventFromTE(TE, 'Port1In', cueLickWindow, TE.Cue);
TE.usLicks = countEventFromTE(TE, 'Port1In', usLickWindow, TE.Us);

%% generate trial lookups for different combinations of conditions
    validTrials = filterTE(TE, 'reject', 0);
    highValueTrials = filterTE(TE, 'trialType', 1:3, 'reject', 0);
    lowValueTrials = filterTE(TE, 'trialType', 4:6, 'reject', 0);
    uncuedTrials = filterTE(TE, 'trialType', 7:9, 'reject', 0);    
    rewardTrials = filterTE(TE, 'trialOutcome', 1, 'reject', 0);
    punishTrials = filterTE(TE, 'trialOutcome', 2, 'reject', 0);    
    omitTrials = filterTE(TE, 'trialOutcome', 3, 'reject', 0);
    trialTypes = 1:9;
    trialsByType = cell(size(trialTypes));
    for counter = 1:length(trialTypes)
        trialsByType{counter} = filterTE(TE, 'trialType', trialTypes(counter), 'reject', 0);
    end
    
    
    
%% compute residuals w.r.t. average and add to TE
% TE.phPeak_cs = bpCalcPeak_dFF(TE.Photometry, 1, [0 2], TE.Cue, 'method', 'mean');
% TE.phPeak_us = bpCalcPeak_dFF(TE.Photometry, 1, [0 0.5], TE.Us, 'method', 'mean');

% make sure cue comes on at 4s and Us at 7s
%     for trial = 1:length(TE.epoch)        
%         cueOn = TE.Cue{trial}(1) - TE.PreCsRecording{trial}(1);
%         UsOn = TE.Us{trial}(1) - TE.PreCsRecording{trial}(1);
%         if cueOn < 3.5 || cueOn > 4.5 || UsOn < 6.5 || UsOn > 7.5
%             disp(['problem with trial ' num2str(trial)]);
%         end
%     end
%     %
    trialTypeAverages = struct('avg', 'avgSEM', 'phPeak_cs', 'phPeak_us');
%     cueTimes = cellfun(@(x,y) x(1) - y(1), TE.Cue, TE.PreCsRecording);
%     usTimes = cellfun(@(x,y) x(1) - y(1), TE.Us, TE.PreCsRecording);   
    % figure out how long the baseline is since this varies across animals
    % (3 or 4)
    bl = round(TE.Cue{1}(1) - TE.PreCsRecording{1}(1));
    cueZero = -bl;
    usZero = -bl - 3;

    cuePoints = [bpX2pnt(cuePhWindow(1), 20, cueZero) bpX2pnt(cuePhWindow(2), 20, cueZero)];
    usPoints = [bpX2pnt(usPhWindow(1), 20, usZero) bpX2pnt(usPhWindow(2), 20, usZero)];
    TE.phResidual_cs = zeros(size(TE.filename));
    TE.phResidual_us = zeros(size(TE.filename));
    for counter = 1:length(trialsByType)
        [trialTypeAverages(counter).avg, trialTypeAverages(counter).avgSEM] = phAverageFromTE(TE, trialsByType{counter}, 1);
        trialTypeAverages(counter).phPeak_cs = nanmean(trialTypeAverages(counter).avg(cuePoints(1):cuePoints(2)));
        trialTypeAverages(counter).phPeak_us = nanmean(trialTypeAverages(counter).avg(usPoints(1):usPoints(2)));
        for trial = find(trialsByType{counter})
            TE.phResidual_cs(trial) = TE.phPeak_cs.data(trial) - trialTypeAverages(counter).phPeak_cs;
            TE.phResidual_us(trial) = TE.phPeak_us.data(trial) - trialTypeAverages(counter).phPeak_us;            
        end
    end
    %%
    
    ensureFigure('cuedOutcome_correlation_plots', 1);
    
    % us vs cue licks
    subplot(2,2,1);
    cuedRewardTrials = find(trialsByType{1} | trialsByType{4});  
    cuedRewardTrials = setdiff(cuedRewardTrials, cuedRewardTrials(find(isnan(TE.csLicks.count(cuedRewardTrials))))); % get rid of NaNs in lick counts
    hiRewardTrials = find(trialsByType{1});  
    hiRewardTrials = setdiff(hiRewardTrials, hiRewardTrials(find(isnan(TE.csLicks.count(hiRewardTrials)))));
    lowRewardTrials = find(trialsByType{4});  
    lowRewardTrials = setdiff(lowRewardTrials, lowRewardTrials(find(isnan(TE.csLicks.count(lowRewardTrials)))));
    
    scatter(TE.csLicks.count(cuedRewardTrials), TE.phResidual_us(cuedRewardTrials), 'b', '.'); hold on;
    fo = fitoptions('poly1', 'Exclude', TE.csLicks.count(cuedRewardTrials) > 50);%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);

    fob = fit(TE.csLicks.count(cuedRewardTrials), TE.phResidual_us(cuedRewardTrials), 'poly1', fo); 
    plot(fob,'predfunc'); legend off;
    
    set(gca, 'FontSize', 10); xlabel('Cue licks'); ylabel('Rew. (dFF residual)');
    title('Cued reward, hiVal, loVal, combined', 'FontSize', 10); 
    set(gca, 'XLim', [0 20]);
    %     set(gca, 'YLim', [-0.005, 0.01]);
    
    % us vs cue photometry
    subplot(2,2,2);
    scatter(TE.phResidual_cs(cuedRewardTrials), TE.phResidual_us(cuedRewardTrials), 'b', '.'); hold on;
    fo = fitoptions('poly1');%, 'Upper', [0, Inf], 'Lower', [-Inf, 0]);


    fob = fit(TE.phResidual_cs(cuedRewardTrials), TE.phResidual_us(cuedRewardTrials), 'poly1', fo); 
    plot(fob,'predfunc'); legend off;
    textBox(TE.filename{1}(1:7), gca, [0.2 0.95]); 
    set(gca, 'FontSize', 10); xlabel('Cue (dFF residual)'); ylabel('Reward (dFF residual)');
    title('Cued reward, hiVal, loVal, combined', 'FontSize', 10); 
%     set(gca, 'YLim', [-0.005, 0.01]);
    
    % us vs cue licks, high and low separate    
    subplot(2,2,3);
    scatter(TE.csLicks.count(hiRewardTrials), TE.phPeak_us.data(hiRewardTrials), 'k', '.'); hold on;
    scatter(TE.csLicks.count(lowRewardTrials), TE.phPeak_us.data(lowRewardTrials), 'r', '.');    
    set(gca, 'FontSize', 10); xlabel('Cue licks'); ylabel('Reward (dFF)');
    title('Cued reward', 'FontSize', 10);
%     legend({'hival', 'loval'}, 'Location', 'northwest', 'FontSize', 10); legend('boxoff');
%     set(gca, 'YLim', [-0.005, 0.015]);
    set(gca, 'XLim', [0 20]);

    % us vs cue, high and low separate   
    subplot(2,2,4);
    scatter(TE.phPeak_cs.data(hiRewardTrials), TE.phPeak_us.data(hiRewardTrials), 'k', '.'); hold on;
    scatter(TE.phPeak_cs.data(lowRewardTrials), TE.phPeak_us.data(lowRewardTrials), 'r', '.');    
    set(gca, 'FontSize', 10); xlabel('Cue (dFF)'); ylabel('Reward (dFF)');
    title('Cued reward', 'FontSize', 10);
%     legend({'hival', 'loval'}, 'Location', 'northwest', 'FontSize', 10); legend('boxoff');
%     set(gca, 'YLim', [-0.005, 0.015]);
if saveOn    
    saveas(gcf, fullfile(savepath, 'cuedOutcome_correlation_plots.fig'));
    saveas(gcf, fullfile(savepath, 'cuedOutcome_correlation_plots.jpg'));
end
    
if saveOn
    save(fullfile(savepath, 'TE.mat'), 'TE');
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
end    
    