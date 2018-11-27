function DB = dbRegisterAnimal(DB, animal, delete)
% ensure that animal folders are created in database
% save and return the updated DB file for the experiment

if nargin < 3
    delete = 0;
end

if ~delete
    ensureDirectory(fullfile(DB.path, 'animals', [animal filesep]));
end

if ~iscell(DB.animals)
    DB.animals = {};
end


matches = strcmp(DB.animals, animal);
if ~delete
    if ~any(matches)
        DB.animals{end+1} = animal;
    end
else
    DB.animals = DB.animals(~matches);
end

if ~delete
    fprintf('*** animal %s registered in experiment %s\n', animal, DB.name);
else
    fprintf('*** animal %s removed from experiment %s\n', animal, DB.name);    
end
save(fullfile(DB.path, 'DB.mat'), 'DB');