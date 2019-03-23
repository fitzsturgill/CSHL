
basepath = 'Z:\SummaryAnalyses\ChAT_PE_Manuscript';
savePath = fullfile(basepath, 'pie_charts');
ensureDirectory(savePath);



saveName = 'pie_surprise_modulation';
ensureFigure(saveName, 1);

cue_mod = [.001 8];
reward = [.001 8];
surprise_mod = [6 2];
subplot(3,1,1);
pie(reward);
title('reward');
subplot(3,1,2);
p = pie(cue_mod);
title('cue');
subplot(3,1,3);
pie(surprise_mod);
title('RPE');

formatFigurePublish('size', [1 3]);


if saveOn 
    export_fig(fullfile(savePath, saveName), '-eps');
end
