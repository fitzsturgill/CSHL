% cuedOutcome_Odor_complete summary analysis script 9/17
% includes separate analysis for phasic and sustained components to odor
% cue response

files = {...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_34', 'ccomplete_outcomeSummary_ChAT_34.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_35', 'ccomplete_outcomeSummary_ChAT_35.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_39', 'ccomplete_outcomeSummary_ChAT_39.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_42', 'ccomplete_outcomeSummary_ChAT_42.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_26', 'ccomplete_outcomeSummary_ChAT_26.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_32', 'ccomplete_outcomeSummary_ChAT_32.mat';...
    'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\ChAT_37', 'ccomplete_outcomeSummary_ChAT_37.mat'}; % this one last because of weird emergence and dissapearence of reward responses

for counter = 1:size(files, 1)
    fileName = files{counter, 2}(1:end-4);
    temp = load(fullfile(files{counter, 1}, files{counter, 2}));
    if counter == 1
        sumOutcomeData = temp.ccomplete_outcomeSummary; % copy it first time, then add
        sumOutcomeData(1).filename = {fileName};
        fnames = fieldnames(temp.ccomplete_outcomeSummary)';     % field names must all be the same, (avg, n, std, etc.)       
    else
        thisData = temp.ccomplete_outcomeSummary;
        sumOutcomeData(1).filename{end + 1} = fileName;
        for type = 1:9
            for f = fnames
                sumOutcomeData(type).(f{:})(end + 1) = thisData(type).(f{:});
            end
        end          
    end
end

%% reward bar graph
savepath = 'C:\Users\Adam\Dropbox\KepecsLab\_Fitz\CCN\CCN_Talk\';
saveOn = 1;
ensureFigure('cuedOutcome_reward_barGraph', 1);
offset = 0.1;
bar(1, mean(sumOutcomeData(1).avg), 'b'); hold on; 
bar(2, mean(sumOutcomeData(4).avg), 'r'); 
bar(3, mean(sumOutcomeData(7).avg), 'g'); 
errorbar(1, mean(sumOutcomeData(1).avg), mean(sumOutcomeData(1).avg)/sqrt(length(sumOutcomeData(1).avg)), 'b.', 'linewidth', 2);
errorbar(2, mean(sumOutcomeData(4).avg), mean(sumOutcomeData(4).avg)/sqrt(length(sumOutcomeData(4).avg)), 'r.', 'linewidth', 2);
errorbar(3, mean(sumOutcomeData(7).avg), mean(sumOutcomeData(7).avg)/sqrt(length(sumOutcomeData(7).avg)), 'g.', 'linewidth', 2);
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'High V.', 'Low V.', 'Uncued'}, 'YLim', [0 1.2], 'XLim', [0.5 3.5]);
ylabel('\Delta ZScore');
formatFigureTalk([4 2]);   
if saveOn    
    saveas(gcf, fullfile(savepath, 'cuedOutcome_reward_barGraph.fig'));
    saveas(gcf, fullfile(savepath, 'cuedOutcome_reward_barGraph.jpg'));    
    saveas(gcf, fullfile(savepath, 'cuedOutcome_reward_barGraph.meta'));                   
end

%% punish bar graph
ensureFigure('cuedOutcome_punish_barGraph', 1);
offset = 0.1;
bar([1 2 3], [mean(sumOutcomeData(2).avg) mean(sumOutcomeData(5).avg) mean(sumOutcomeData(8).avg)], 'r'); hold on;
errorbar([1 2 3], [mean(sumOutcomeData(2).avg) mean(sumOutcomeData(5).avg) mean(sumOutcomeData(8).avg)],...
    [mean(sumOutcomeData(2).avg)/sqrt(length(sumOutcomeData(2).avg)),...
    mean(sumOutcomeData(5).avg)/sqrt(length(sumOutcomeData(5).avg)),...
    mean(sumOutcomeData(8).avg)/sqrt(length(sumOutcomeData(8).avg))], 'r.', 'linewidth', 2);
set(gca, 'XTickLabel', {'Uncued', 'Low V.', 'High V.'});

