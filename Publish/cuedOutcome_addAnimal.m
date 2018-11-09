function cuedOutcome_addAnimal(animal)% cuedOutcome_addAnimal
DB = dbLoadExperiment('cuedOutcome');
sourcePath = 'Z:\SummaryAnalyses\CuedOutcome_Odor_Complete\';

DB = dbRegisterAnimal(DB, animal);
savepath = dbGetAnimalPath(DB, animal);

success = copyfile(fullfile([sourcePath animal filesep], 'TE.mat'), fullfile(savepath, 'TE.mat'));
if success
    disp(['*** Saved: ' fullfile(savepath, 'TE.mat')]);
else
    warning('*** adding %s to cuedOucome database failed ***', animal);
end