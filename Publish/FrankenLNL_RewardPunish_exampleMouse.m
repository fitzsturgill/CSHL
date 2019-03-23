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


%% example traces
% rewarding subset
window = [-4 7];
showTheseR = find(Odor2Valve1Trials & rewardTrials & TE.sessionIndex == sessionIndex);
[~, rixR] = sort(rand(size(showTheseR)));



saveName = 'PE_BLA_exampleMouse_traces';  
fig = ensureFigure(saveName, 1);    


ax = subplot(1,1,1);
plot(TE.Photometry.xData, TE.Photometry.data(1).dFF(showTheseR(rixR(1:10)), :)', 'k', 'LineWidth', 0.15); set(gca, 'XLim', window);
addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [2.9 3.1]);
% set(gca, 'YTickLabel', {}); set(gca, 'XTickLabel', {}); 
set(gca, 'Visible', 'off');
% title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Left']);


% 
% ax(2) = subplot(1,2,2);
% plot(TE.Photometry.xData, TE.Photometry.data(2).dFF(showTheseR(rixR(1:10)), :)', 'k', 'LineWidth', 0.15); set(gca, 'XLim', window);
% addStimulusPatch(gca, [0 1]); addStimulusPatch(gca, [2.9 3.1]);
% % title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Left']);
% % set(gca, 'YTickLabel', {}); set(gca, 'XTickLabel', {}); 
% set(gca, 'Visible', 'off');


% aversive subset
showTheseA = find(Odor2Valve2Trials & punishTrials & TE.sessionIndex == sessionIndex);
[~, rixA] = sort(rand(size(showTheseA)));

formatFigurePublish('size', [1.6 1.1]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end

        
%% averages

fdField = 'ZS';
saveName = 'PE_BLA_exampleMouse_avgs';  
h=ensureFigure(saveName, 1); 


linecolors = [0 0 1; 0 0 0; 0 1 1];            

subplot(2, 1, 1);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Left']);
addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [2.9 3.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {trialsByType{1} & TE.sessionIndex == sessionIndex, trialsByType{2} & TE.sessionIndex == sessionIndex, trialsByType{5} & TE.sessionIndex == sessionIndex,}, 1,...
    'FluorDataField', fdField, 'window', [-4, 7], 'cmap', linecolors); %high value, reward

% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northwest'); legend('boxoff');
ylabel('(\fontsize{12}\sigma\fontsize{8}-baseline)');  set(gca, 'XLim', [-4 7]);

subplot(2, 1, 2);
set(gca, 'YLim', [-2 10]);
tcolor = mycolors('chat');
title(['\color[rgb]{' sprintf('%.4f,%.4f,%.4f', tcolor(1), tcolor(2), tcolor(3)) '}Right']);
addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7], 0.4);  addStimulusPatch(gca, [2.9 3.1], '', [0.7 0.7 0.7], 0.4);
[ha, hl] = phPlotAverageFromTE(TE, {trialsByType{1} & TE.sessionIndex == sessionIndex, trialsByType{2} & TE.sessionIndex == sessionIndex, trialsByType{5} & TE.sessionIndex == sessionIndex,}, 2,...
    'FluorDataField', fdField, 'window', [-4, 7], 'cmap', linecolors); %high value, reward

% legend(hl, {'omit', 'cued', 'uncued'}, 'Location', 'northwest'); legend('boxoff');
ylabel('\fontsize{8}Fluor. '); xlabel('time from cue (s)'); set(gca, 'XLim', [-4 7]);


formatFigurePublish('size', [1.5 2]);

if saveOn 
    export_fig(fullfile(savepath, saveName), '-eps');
end
