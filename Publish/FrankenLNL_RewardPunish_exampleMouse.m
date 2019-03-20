% FrankenLNL_RewardPunish_exampleMouse

DB = dbLoadExperiment('FrankenLNL_RewardPunish');
savepath = fullfile(DB.path, ['pooled' filesep 'figure']);
ensureDirectory(savepath);
smoothWindow = 1;
saveOn = 1;
animal = 'ACh_7';

photometryField = 'Photometry';
fdField = 'ZS';

success = dbLoadAnimal(DB, animal); % load TE and trial lookups

%% % find the session index for example session with air puff
sessionName = 'ACh_7_FrankenLNL_4odors_Feb27_2019_Session1.mat';
matches = strcmp(TE.filename, sessionName);
sessionIndex = unique(TE.sessionIndex(matches));

        
%% averages

fdField = 'ZS';
saveName = sprintf('%s_phAvgs_%s', animal, fdField);  
h=ensureFigure(saveName, 1); 


linecolors = [1 0 0; 0 0 1; 0 1 1; 0 1 0];            

subplot(1, 2, 1);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Ach.']);
addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [2.9 3.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {neutralTrials & csPlusTrials & hitTrials, rewardTrials & csPlusTrials & hitTrials, csMinusTrials & CRTrials & rewardTrials, uncuedReward}, 1,...
    'FluorDataField', fdField, 'window', [-4, 7], 'cmap', linecolors); %high value, reward

% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northwest'); legend('boxoff');
ylabel('\fontsize{8}Fluor. (\fontsize{12}\sigma\fontsize{8}-baseline)'); xlabel('time from cue (s)'); set(gca, 'XLim', [-4 7]);



subplot(1, 2, 2);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('dat');
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Dop.']);
addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [2.9 3.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {neutralTrials & csPlusTrials & hitTrials, rewardTrials & csPlusTrials & hitTrials, csMinusTrials & CRTrials & rewardTrials, uncuedReward}, 2,...
    'FluorDataField', fdField, 'window', [-4, 7], 'cmap', linecolors); %high value, reward
% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northoutside'); legend('boxoff');
ylabel('\fontsize{8}Fluor. (\fontsize{12}\sigma\fontsize{8}-baseline)'); xlabel('time from cue (s)'); set(gca, 'XLim', [-4 7]);

formatFigurePublish('size', [3 1.1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end
