%% grand averages

files = {...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_34', 'averages_ChAT_34.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_35', 'averages_ChAT_35.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_39', 'averages_ChAT_39.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_42', 'averages_ChAT_42.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_26', 'averages_ChAT_26.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_32', 'averages_ChAT_32.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_37', 'averages_ChAT_37.mat'}; % this one last because of weird emergence and dissapearence of reward responses

dataFields = {'cue', 'full'}; subFields = {'photometry', 'licks'}; subSubFields = {'avg', 'xData'};
for counter = 1:size(files, 1)
    fileName = files{counter, 2}(1:end-4);
    load(fullfile(files{counter, 1}, files{counter, 2}));    % load avgData into workspace
    if counter == 1
        avgData = data;
    else
        for dfi = 1:length(dataFields)
            for sdfi = 1:length(subFields)
                for ssdfi = 1:length(subSubFields)
                    avgData.(dataFields{dfi}).(subFields{sdfi}).(subSubFields{ssdfi})(counter, :, :) = data.(dataFields{dfi}).(subFields{sdfi}).(subSubFields{ssdfi});
                end
            end
        end
    end
end

%% plot grand average cue lick and photometry responses, normalized by high value cue area
% 3rd dimension - first chunk = high value, second = low value, third =
% uncued
cueLicks_norm = avgData.cue.licks.avg;
% baseline licks
cuePh_norm = avgData.cue.licks.avg;
nSessions = size(cueLicks_norm, 1);
xData_licks = squeeze(avgData.cue.licks.xData(1,:,1)); zp_licks = find(xData_licks > 0, 1);
xData_ph = squeeze(avgData.cue.photometry.xData(1,:,1)); zp_ph = find(xData_ph > 0, 1);
% baseline licks
cueLicks_baselined = bsxfun(@minus, cueLicks_norm, mean(cueLicks_norm(:,1:zp_licks - 1, :), 2));
lickD = zeros(nSessions, 1);
phD = zeros(nSessions, 1);
for counter = 1:nSessions
    lickD = trapz(xData_licks(zp_licks:end), squeeze(cueLicks_baselined(counter, zp_licks:end, 1))); % lick denominator for normalization
    phD = trapz(xData_ph(zp_ph:end), squeeze(avgData.cue.photometry.avg(counter, zp_ph:end, 1))); % photometry denominator for normalization
end
cueLicks_norm = bsxfun(@rdivide, cueLicks_baselined, lickD);
cuePh_norm = bsxfun(@rdivide, avgData.cue.photometry.avg, phD);

%%
figName = 'Cue_grandAverage';
ensureFigure(figName, 1);
% mean
subplot(1,2,1);
hl = plot(xData_licks, squeeze(mean(cueLicks_norm))', '--'); hold on;
hl(1).Color = 'b'; hl(2).Color = 'r'; hl(3).Color = 'g';
hl = plot(xData_ph, squeeze(mean(cuePh_norm))', '-'); hold on;
hl(1).Color = 'b'; hl(2).Color = 'r'; hl(3).Color = 'g';
set(gca, 'XLim', [-2 3]); xlabel('time from cue (s)'); ylabel('normalized by high value integral (0-3s)'); title('mean');
legend({'Lick, high', 'Lick, low', 'Lick, uncued', 'Ph, high', 'Ph, low', 'Ph, uncued'}, 'Location', 'northwest', 'Box', 'off');
% median
subplot(1,2,2);
hl = plot(xData_licks, squeeze(median(cueLicks_norm))', '--'); hold on;
hl(1).Color = 'b'; hl(2).Color = 'r'; hl(3).Color = 'g';
hl = plot(xData_ph, squeeze(median(cuePh_norm))', '-'); hold on;
hl(1).Color = 'b'; hl(2).Color = 'r'; hl(3).Color = 'g';
set(gca, 'XLim', [-2 3]);xlabel('time from cue (s)'); ylabel('normalized by high value integral (0-3s)'); title('median');
legend({'Lick, high', 'Lick, low', 'Lick, uncued', 'Ph, high', 'Ph, low', 'Ph, uncued'}, 'Location', 'northwest', 'Box', 'off');
formatFigureEvernote([10 5]);

%% bounded line version
savepath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\grandAverages\';
saveOn = 1;
figName = 'Cue_grandAverage_SEM';
ensureFigure(figName, 1); axes; set(gca, 'YLim', [-0.2 1.2]);
addStimulusPatch(gca, [0 1], '', [0.7 0.7 0.7]);
% mean
cmap = [0 0 1; 1 0 0; 0 1 0];
cmap_licks = cmap * 0.33;
cmap_ph = cmap;
cueLicks_norm_avg = squeeze(mean(cueLicks_norm));
cueLicks_norm_sem = squeeze(std(cueLicks_norm)) / sqrt(size(cueLicks_norm, 1));

[lines, patches] = boundedline(xData_licks', cueLicks_norm_avg, permute(cueLicks_norm_sem, [1 3 2]), 'cmap', cmap_licks); hold on;

cuePh_norm_avg = squeeze(mean(cuePh_norm));
cuePh_norm_sem = squeeze(std(cuePh_norm)) / sqrt(size(cuePh_norm, 1));
[hl, hp] = boundedline(xData_ph', cuePh_norm_avg, permute(cuePh_norm_sem, [1 3 2]), 'cmap', cmap_ph); 
patches = [patches; hp];
set(gca, 'XLim', [-2 3]); xlabel('time from cue (s)'); ylabel('normalized by high value integral (0-3s)');
legend(patches, {'Lick, high', 'Lick, low', 'Lick, uncued', 'Ph, high', 'Ph, low', 'Ph, uncued'}, 'Location', 'northwest', 'Box', 'off');

formatFigureEvernote([7 5]);
if saveOn    
    saveas(gcf, fullfile(savepath, [figName '.fig']));
    saveas(gcf, fullfile(savepath, [figName '.jpg']));
    saveas(gcf, fullfile(savepath, [figName '.meta']));                 
end

%% plot grand average full lick and photometry responses, normalized by high value cue 90%
% 3rd dimension - first chunk = high value, second = low value, third =
% uncued
fullLicks_norm = avgData.full.licks.avg;
% baseline licks
fullPh_norm = avgData.full.licks.avg;
nSessions = size(fullLicks_norm, 1);
xData_licks = squeeze(avgData.full.licks.xData(1,:,1)); zp_licks = find(xData_licks > 0, 1);
xData_ph = squeeze(avgData.full.photometry.xData(1,:,1)); zp_ph = find(xData_ph > 0, 1);
% baseline licks
fullLicks_baselined = bsxfun(@minus, fullLicks_norm, mean(fullLicks_norm(:,1:zp_licks - 1, :), 2));
lickD = zeros(nSessions, 1);
phD = zeros(nSessions, 1);
for counter = 1:nSessions
    lickD = trapz(xData_licks(zp_licks:end), squeeze(fullLicks_baselined(counter, zp_licks:end, 1))); % lick denominator for normalization
    phD = trapz(xData_ph(zp_ph:end), squeeze(avgData.full.photometry.avg(counter, zp_ph:end, 1))); % photometry denominator for normalization
end
fullLicks_norm = bsxfun(@rdivide, fullLicks_baselined, lickD);
fullPh_norm = bsxfun(@rdivide, avgData.full.photometry.avg, phD);




        
        


