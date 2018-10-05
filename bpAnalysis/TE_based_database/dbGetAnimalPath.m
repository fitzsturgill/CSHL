function path = dbGetAnimalPath(DB, animal)

    % return animal path for experiment (as defined by the experiment DB
    % variable)
    
    path = fullfile(DB.path, 'animals', [animal filesep]);
    
    if exist(path) == 7
        return
    else
        path = '';
        return
    end
    