% CuedOutcome_QPE_examples
DB = dbLoadExperiment('cuedOutcome');
savepath = fullfile(DB.path, 'figure');
ensureDirectory(savepath);
saveOn = 1;


% photometryField = 'Photometry';
% fdField = 'ZS';


%% Lick histogram for graded value task
figSize = [1.7 0.9];
animal = 'ChAT_39';
xlim = [-3 3];
success = dbLoadAnimal(DB, animal); % load TE and trial lookups
 % cue types
    varargin = {'trialNumbering', 'consecutive',...
        'window', [-6 4], 'zeroField', 'Us', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'m', 'b', 'k'}};
    saveName = 'Lick_Hist_CCN';
    ensureFigure(saveName, 1); axes;

    [ha, hl] = plotEventAverageFromTE(TE, trialsByType([1 4 7]), 'Port1In', varargin{:});
    set(gca, 'XLim', xlim, 'YLim', [-2 15], 'XTick', [-3 0 3]);
    addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7]) 
    addStimulusPatch(gca, [-0.1 0.1], '');
%     legend(hl, {'\color{blue} high value', '\color{red} low value', '\color{green} uncued'}, 'FontSize', 6, 'Interpreter', 'tex', 'Position', [0.25 0.6 0.3 0.4]); legend('boxoff');    
%     ylabel('Licks (Hz)'); xlabel('Time from reinforcement (s)');
    set(gca, 'YTick', [0 10]);
    formatFigurePublish('size', figSize);
    if saveOn 
        print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
        export_fig(fullfile(savepath, saveName), '-eps');
    end    
    
    
%% avgs
figSize = [1.7 0.9];


animal = 'ChAT_42';
xlim = [-4 4];
success = dbLoadAnimal(DB, animal); % load TE and trial lookups

% reward and punish avgs
saveName = 'CuedOutcome_avgs_complete';
ensureFigure(saveName, 1); axes;
[ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
    'window', [xlim(1) 0], 'linespec', {'b', 'm'}, 'FluorDataField', 'ZS'); hold on;

[ha, hl] = phPlotAverageFromTE(TE, {rewardTrials, punishTrials}, 1,...
    'window', [0 xlim(2)], 'cmap', [mycolors('reward'); mycolors('puff')], 'FluorDataField', 'ZS');
hl = [hla hl]; 
addStimulusPatch(gca, [-3 -2], '', [0.8 0.8 0.8], 0.5);
addStimulusPatch(gca, [-0.1 0.1], '', [0.8 0.8 0.8], 0.5);
formatFigurePublish('size', figSize);
set(gca, 'XLim', xlim);
if saveOn
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end  
    
% reward avgs

saveName = 'CuedOutcome_avgs_reward';
ensureFigure(saveName, 1); axes;
[ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials & rewardTrials, highValueTrials & rewardTrials}, 1,...
    'window', xlim, 'linespec', {'b', 'm'}, 'FluorDataField', 'ZS');
addStimulusPatch(gca, [-3 -2], '', [0.8 0.8 0.8], 0.5);
addStimulusPatch(gca, [-0.1 0.1], '', [0.8 0.8 0.8], 0.5);
formatFigurePublish('size', figSize);
set(gca, 'XLim', xlim);
if saveOn
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end  
% punish avgs

saveName = 'CuedOutcome_avgs_punish';
ensureFigure(saveName, 1); axes;
[ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials & punishTrials, highValueTrials & punishTrials}, 1,...
    'window', xlim, 'linespec', {'b', 'm'}, 'FluorDataField', 'ZS');
addStimulusPatch(gca, [-3 -2], '', [0.8 0.8 0.8], 0.5);
addStimulusPatch(gca, [-0.1 0.1], '', [0.8 0.8 0.8], 0.5);
formatFigurePublish('size', figSize);
set(gca, 'XLim', xlim);
if saveOn
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end
%% lets try cue, reward, punish separately

figSize = [1.7 0.9];


animal = 'ChAT_42';

success = dbLoadAnimal(DB, animal); % load TE and trial lookups

% cue only
window = [-1 3];
saveName = 'CuedOutcome_avgs_cueOnly';
ensureFigure(saveName, 1); axes;
[ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials, highValueTrials}, 1,...
    'window', window, 'zeroTimes', TE.Cue, 'linespec', {'b', 'm'}, 'FluorDataField', 'ZS'); hold on;


addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);

formatFigurePublish('size', figSize);
set(gca, 'XLim', window);
if saveOn
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end  

% reward only
window = [0 4];
saveName = 'CuedOutcome_avgs_rewardOnly';
ensureFigure(saveName, 1); axes;
[ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials & rewardTrials, highValueTrials & rewardTrials}, 1,...
    'window', window, 'zeroTimes', TE.Us, 'linespec', {'b', 'm'}, 'FluorDataField', 'ZS'); hold on;


addStimulusPatch(gca, [0 0.1], '', [0.8 0.8 0.8], 0.5);

formatFigurePublish('size', figSize);
set(gca, 'XLim', window);
if saveOn
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end  

% punish only
window = [0 4];
saveName = 'CuedOutcome_avgs_punishOnly';
ensureFigure(saveName, 1); axes;
[ha, hla] = phPlotAverageFromTE(TE, {lowValueTrials & punishTrials, highValueTrials & punishTrials}, 1,...
    'window', window, 'zeroTimes', TE.Us, 'linespec', {'b', 'm'}, 'FluorDataField', 'ZS'); hold on;


addStimulusPatch(gca, [0 0.1], '', [0.8 0.8 0.8], 0.5);

formatFigurePublish('size', figSize);
set(gca, 'XLim', window);
if saveOn
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
end  
%% combined licking and photometry averages from ChAT_42


    load('Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_42\TE.mat');
    cuedOutcome_Conditions;
    saveName = 'CuedOutcome_example_averages_combined_ChAT_39';
    window = [-1.5 5];
    ensureFigure(saveName, 1); 
        varargin = {'window', [window(1) 3], 'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording',...
        'linespec', {'m', 'b'}};
    lickAvg = eventAverageFromTE(TE, {highValueTrials & rewardTrials, lowValueTrials & rewardTrials}, 'Port1In', varargin{:});
    ax = axes; hold on; yyaxis right
    ll = plot(lickAvg.xData, lickAvg.Avg(1,:), '--k');
    plot(lickAvg.xData, lickAvg.Avg(1,:), '--m', lickAvg.xData, lickAvg.Avg(2,:), '--b');
    ax.YColor = [0 0 0]; ylabel('Licks/s');
%     set(gca, 'YLim', [-1 7]);
    yyaxis left; ax.YColor = [0 0 0]; 
    [ha, hl] = phPlotAverageFromTE(TE, {highValueTrials & rewardTrials, lowValueTrials & rewardTrials}, 1,...
        'window', window, 'linespec', {'m', 'b'}, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue, 'alpha', 0);
%     set(hl, 'LineWidth', 2);    
    set(gca, 'XLim', window);
    addStimulusPatch(gca, [0 1], '', [0.8 0.8 0.8], 0.5);
    addStimulusPatch(gca, [2.9 3.1], '', [0.8 0.8 0.8], 0.5); 
%     legend([hl ll(1)], {'\color{blue} high value', '\color{red} low value', '\color{green} uncued', ...
%         'licking'}, 'Location', 'northwest', 'FontSize', 7, 'Interpreter', 'tex'); legend('boxoff');
    ylabel('Fluor. (\sigma-bl.)'); xlabel('Time from cue (s)');
    formatFigurePublish('size', [2 1]);    

if saveOn
%     export_fig(fullfile(savepath, saveName), '-transparent', '-pdf');
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
    saveas(gcf, fullfile(savepath, [saveName '.epsc']));   
end    


%% lick-aligned rasters from ChAT_42 version #1

animal = 'ChAT_42';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups

photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
channels = [1];
climfactor = 2;
lickTickLineWidth = 0.3; % 0.3 works well when printed out given current sizing

fontsize = 8;

savepath = fullfile(DB.path, 'figure', filesep);
ensureDirectory(savepath);

xwindow = [-2 5];

saveName = sprintf('%s_lickAligned_rasters_%s_%s_v1', animal, photometryField, fdField);  
ensureFigure(saveName, 1); 

lickOnsets = TE.lickLatency_cs(highValueTrials & rewardTrials);
lickOnsets = sort(lickOnsets);    
subplot(1,4,1);
eventRasterFromTE(TE, highValueTrials & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.lickLatency_cs, 'LineWidth', lickTickLineWidth);
set(gca, 'XLim', xwindow);
% title('high value'); 
ylabel('trial # (sorted)');

subplot(1,4,2);
phRasterFromTE(TE, highValueTrials & rewardTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'sortValues', TE.lickLatency_cs, 'zeroTimes', TE.Cue, 'window', xwindow); % 'CLimFactor', CLimFactor,
line(lickOnsets, (1:sum(highValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1);    
xlabel('Time frome odor (s)');

lickOnsets = TE.lickLatency_cs(lowValueTrials & rewardTrials);
lickOnsets = sort(lickOnsets);        
subplot(1,4,3);
eventRasterFromTE(TE, lowValueTrials & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.lickLatency_cs, 'LineWidth', lickTickLineWidth);
set(gca, 'XLim', xwindow);
% title('low value');

subplot(1,4,4);
phRasterFromTE(TE, lowValueTrials & rewardTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'sortValues', TE.lickLatency_cs, 'zeroTimes', TE.Cue, 'window', xwindow); % 'CLimFactor', CLimFactor,             
line(lickOnsets, (1:sum(lowValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1);
axs = findobj(gcf, 'Type', 'axes');
set(axs, 'FontSize', fontsize);
set(axs([1 3]), 'YTick', []);

formatFigurePublish('size', [4 1.2]);

if saveOn
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%     export_fig(fullfile(savepath, saveName), '-eps');
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
end  

%% lick-aligned rasters from ChAT_42:  VERSION #2
% only show high value
animal = 'ChAT_42';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups

photometryField = 'Photometry';
fdField = 'ZS';
saveOn = 1;
channels = [1];
climfactor = 2;
lickTickLineWidth = 0.1; % 0.3 works well when printed out given current sizing

fontsize = 7;

savepath = fullfile(DB.path, 'figure', filesep);
ensureDirectory(savepath);

xwindow = [-2 5];

saveName = sprintf('%s_lickAligned_rasters_%s_%s', animal, photometryField, fdField);  
ensureFigure(saveName, 1); 

lickOnsets = TE.lickLatency_cs(highValueTrials & rewardTrials);
lickOnsets = sort(lickOnsets);    
subplot(1,2,1);
eventRasterFromTE(TE, highValueTrials & rewardTrials, 'Port1In', 'trialNumbering', 'consecutive',...
    'zeroField', 'Cue', 'startField', 'PreCsRecording', 'endField', 'PostUsRecording', 'sortValues', TE.lickLatency_cs, 'LineWidth', lickTickLineWidth);
set(gca, 'XLim', xwindow);
% title('high value'); 
ylabel('# sorted');
xlabel('Time from odor (s)');

subplot(1,2,2);
phRasterFromTE(TE, highValueTrials & rewardTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'sortValues', TE.lickLatency_cs, 'zeroTimes', TE.Cue, 'window', xwindow); % 'CLimFactor', CLimFactor,
line(lickOnsets, (1:sum(highValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1);    
% colorbar;
xlabel('Time from odor (s)');
set(gca, 'YTick', []);
% lickOnsets = TE.lickLatency_cs(lowValueTrials & rewardTrials);
% lickOnsets = sort(lickOnsets);        
% subplot(1,3,3);
% phRasterFromTE(TE, lowValueTrials & rewardTrials, 1, 'trialNumbering', 'consecutive', 'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', photometryField, 'sortValues', TE.lickLatency_cs, 'zeroTimes', TE.Cue, 'window', xwindow); % 'CLimFactor', CLimFactor,             
% line(lickOnsets, (1:sum(lowValueTrials & rewardTrials))', 'Parent', gca, 'Color', 'r', 'LineWidth', 1);
% set(gca, 'YTick', [0 100 200]);



formatFigurePublish('size', [2.2 0.7]);

if saveOn
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
end  


%% reward and punishment response examples for ChAT_26 (only other possibility being ChAT_39)
animal = 'ChAT_26';
success = dbLoadAnimal(DB, animal); % load TE and trial lookups
figSize = [1 1];
saveName = 'RewardPunish_example_avg';
ensureFigure(saveName, 1);
window = [-1 4];
[ha, hl] = phPlotAverageFromTE(TE, {uncuedTrials & rewardTrials, uncuedTrials & punishTrials}, 1,...
        'window', window, 'linespec', {'b', 'r'}, 'FluorDataField', 'ZS', 'zeroTimes', TE.Us);
xlabel('Time (s)');
ylabel('Fl. (\sigma-bl.)');
set(gca, 'XLim', window);
formatFigurePublish('size', figSize);    

if saveOn
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
    saveas(gcf, fullfile(savepath, [saveName '.fig']));
    saveas(gcf, fullfile(savepath, [saveName '.jpg']));   
end  