function DB = dbLoadExperiment(experiment)
    
    if ~ispref('DB')
        error('database not initialized');
    end
    exps = getpref('DB', 'experiments');
    expIndex = find(strcmp(exps, experiment));
    if expIndex
        paths = getpref('DB', 'paths');        
        load(fullfile(paths{expIndex}, 'DB.mat'));
    else
        error('experiment doesn''t exist');
    end