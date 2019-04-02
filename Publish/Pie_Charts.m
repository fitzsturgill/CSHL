
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


%%

saveName = 'pie_surprise_modulation_GCaMP';
ensureFigure(saveName, 1);

cue_mod = [.001 15];
reward = [.001 15];
punish = [3 12];
surprise_mod = [15 8];
subplot(2,2,1);
pie(reward);
title('reward');
subplot(2,2,2);
pie(punish);
title('punish');
subplot(2,2,3);
p = pie(cue_mod);
title('cue');
subplot(2,2,4);
pie(surprise_mod);
title('RPE');

formatFigurePublish('size', [2 3]);


if saveOn 
    saveas(gcf, fullfile(savePath, [saveName '.fig']));
    saveas(gcf, fullfile(savePath, [saveName '.jpg']));   
    export_fig(fullfile(savePath, saveName), '-eps');
end