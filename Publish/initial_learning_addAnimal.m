function initial_learning_addAnimal(animal)% cuedOutcome_addAnimal
DB = dbLoadExperiment('initial_learning');
sourcePath = 'Z:\SummaryAnalyses\LNL_odor_v2_first2Sessions\';

DB = dbRegisterAnimal(DB, animal);
savepath = dbGetAnimalPath(DB, animal);

success = copyfile(fullfile([sourcePath animal filesep], 'TE.mat'), fullfile(savepath, 'TE.mat'));
if success
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
else
    warning('*** adding %s to cuedOucome database failed ***', animal);
end