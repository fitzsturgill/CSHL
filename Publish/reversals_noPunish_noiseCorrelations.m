% reversals_noPunish_noiseCorrelations

DB = dbLoadExperiment('reversals_noPunish_publish');
savepath = fullfile(DB.path, 'pooled', filesep);
ensureDirectory(savepath);
saveOn = 1;



%% initialize table to hold correlation coefficients for each animal

comp = {'Lick_vs_ACh'; 'Lick_vs_Dop'; 'ACh_vs_Dop'};
acq_p = zeros(3,1);
ext_p = zeros(3,1);
cp_stats = table(comp, acq_p, ext_p);
cp_stats.acq_p(1) = signrank(all_cps(:,1), all_cps(:,2));
cp_stats.acq_p(2) = signrank(all_cps(:,1), all_cps(:,3));
cp_stats.acq_p(3) = signrank(all_cps(:,2), all_cps(:,3));
cp_stats.ext_p(1) = signrank(all_cps(:,4), all_cps(:,5));
cp_stats.ext_p(2) = signrank(all_cps(:,4), all_cps(:,6));
cp_stats.ext_p(3) = signrank(all_cps(:,5), all_cps(:,6));






