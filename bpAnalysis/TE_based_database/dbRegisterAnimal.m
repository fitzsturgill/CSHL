function DB = dbRegisterAnimal(DB, animal)
% ensure that animal folders are created in database
% save and return the updated DB file for the experiment
ensureDirectory(fullfile(DB.path, 'animals', [animal filesep]));

if ~iscell(DB.animals)
    DB.animals = {};
end

if ~any(strcmp(DB.animals, animal))
    DB.animals{end+1} = animal;
end

fprintf('*** animal %s registered in experiment %s\n', animal, DB.name);
save(fullfile(DB.path, 'DB.mat'), 'DB');