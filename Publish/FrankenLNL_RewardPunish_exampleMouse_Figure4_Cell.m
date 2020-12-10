% FrankenLNL_RewardPunish_exampleMouse

DB = dbLoadExperiment('FrankenLNL_RewardPunish');
animal = 'ACh_15';
ch = 2;
savepath = fullfile(DB.path, ['figure' filesep animal]);
ensureDirectory(savepath);
smoothWindow = 1;
saveOn = 1;
window = [-3 3];

photometryField = 'Photometry';
fdField = 'ZS';

success = dbLoadAnimal(DB, animal); % load TE and trial lookups
minRewardLickRate = 2; % at least n Hz licking (during us window for reward trials)
minCueLickRate = 1; % at least n Hz licking (during cs window for reward cued reward trials)

%%  assemble trial sets
    trialSets = struct(...
        'cuedReward', [],...
        'uncuedReward', [],...
        'omitReward', [],...
        'cuedPuff', [],...
        'uncuedPuff', [],...
        'omitPuff', [],...
        'cuedShock', [],...
        'uncuedShock', [],...
        'omitShock', []...
        );
    
    % cuedReward
    trialSets.cuedReward = rewardTrials & Odor2Valve1Trials & (TE.licks_us.rate > minRewardLickRate) & (TE.licks_cs.rate > minCueLickRate) & ismember(TE.BlockNumber, [2 3]); % exclude shock or extinction days
    % uncuedReward
    trialSets.uncuedReward = uncuedReward & (TE.licks_us.rate > minRewardLickRate) & ismember(TE.BlockNumber, [2 3]);       
    % omitReward    
    trialSets.omitReward = neutralTrials & Odor2Valve1Trials & ismember(TE.BlockNumber, [2 3]);           
    % cuedPuff
    trialSets.cuedPuff = punishTrials & Odor2Valve2Trials;
    % uncuedPuff
    trialSets.uncuedPuff = uncuedPunish;
    % omitPuff    
    trialSets.omitPuff = neutralTrials & Odor2Valve2Trials & (TE.BlockNumber == 3);
    % cuedShock    
    trialSets.cuedShock = shockTrials & Odor2Valve2Trials;
    % uncuedShock
    trialSets.uncuedShock = uncuedShock;     
    % omitShock
    trialSets.omitShock = neutralTrials & Odor2Valve2Trials & (TE.BlockNumber == 4);


%% example traces

figSize = [1.6 0.6];
uncuedShockTrials = find(trialSets.uncuedShock);
shockSubset = uncuedShockTrials([1 2 5 6 7 9 11 12 13 15]);

% ensureFigure('test', 1);
% ha=[];
% for counter = 1:17
%     ha(end + 1) = subplot(5,4,counter);
%     plot(TE.PhotometryHF.xData - 1, TE.PhotometryHF.data(1).ZS(uncuedShockTrials(counter), :)');
%     addStimulusPatch(gca, [-2 -1]); addStimulusPatch(gca, [-0.1 0.1]);
%     textBox(num2str(counter));
%     set(gca, 'XLim', [-3 3]);
% end
% sameYScale(ha);
% 
saveName = ['PE_BLA_exampleMouse_traces_shockUncued'];  
fig = ensureFigure(saveName, 1);

ax = subplot(1,1,1);
plot(TE.PhotometryHF.xData - 1, TE.PhotometryHF.data(1).ZS(shockSubset, :)', 'k', 'LineWidth', 0.15); set(gca, 'XLim', window); % xscale shifted relative to reinforcement
addStimulusPatch(gca, [-2 -1]); addStimulusPatch(gca, [-0.1 0.1]);
set(gca, 'YTickLabel', {}); 
set(gca, 'XTickLabel', {}); 
set(gca, 'YTick', [-5 0 5 10], 'XTickLabel', [], 'XLim', [-3 3], 'YLim', [-5 10], 'XTick', [-3 0 3]);
% set(gca, 'YLim', [-0.1 0.3]);

formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%     export_fig(fullfile(savepath, saveName), '-eps');
end


%% rasters shock
figSize = [1.6 0.81];
saveName = 'PE_BLA_exampleMouse_Shock_Rasters_cued';  
ensureFigure(saveName, 1);
climfactor = 10;
clim = [min(TE.PhotometryHF.data(ch).ZS(:)) max(TE.Photometry.data(ch).ZS(:))];
subplot(1,1,1); phRasterFromTE(TE, trialSets.cuedShock, ch, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'PhotometryHF', 'zeroTimes', TE.Us, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
set(gca, 'YTick', [1 20], 'XLim', window, 'YLim', [1 20], 'XTick', [-3 0 3]);
formatFigurePublish('size', figSize);
if saveOn  
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%     export_fig(fullfile(savepath, saveName), '-eps');
end        

saveName = 'PE_BLA_exampleMouse_Shock_Rasters_uncued';  
ensureFigure(saveName, 1);
subplot(1,1,1); phRasterFromTE(TE, trialSets.uncuedShock, ch, 'trialNumbering', 'consecutive',...
    'CLimFactor', climfactor, 'FluorDataField', fdField, 'PhotometryField', 'PhotometryHF', 'zeroTimes', TE.Us, 'window', window, 'showSessionBreaks', 0); % 'CLimFactor', CLimFactor,
set(gca, 'YTick', [1 10], 'XLim', window, 'XTick', [-3 0 3]);
% colorbar;
formatFigurePublish('size', figSize);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%     export_fig(fullfile(savepath, saveName), '-eps');
end        


%% averages, aversive, just do shock
fdField = 'ZS';
saveName = 'PE_BLA_exampleMouse_shock_avgs';  
h=ensureFigure(saveName, 1); 

linecolors = [mycolors('shock_cued'); mycolors('shock'); 0 0 0];            

axes;
% set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
addStimulusPatch(gca, [-3 -2], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [-0.1 0.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {trialSets.cuedShock, trialSets.uncuedShock, trialSets.omitShock}, ch,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window, 'cmap', linecolors); hold on; %high value, reward
% add uncued shock
[ha, hl] = phPlotAverageFromTE(TE, uncuedShock, 1,...
    'zeroTimes', TE.Us, 'FluorDataField', fdField, 'window', window, 'cmap', mycolors('shock')); 
set(gca, 'XLim', window, 'XTick', [-3 0 3]);


formatFigurePublish('size', [2 1]);

if saveOn 
    print(gcf, '-dpdf', fullfile(savepath, [saveName '.pdf']));
%     export_fig(fullfile(savepath, saveName), '-eps');
end

