function dbRegisterExperiment(experiment, path)
    % registers an experiment in preferences and stores path to DB file
    
    % DB preference contains with fields containing cell arrays of size
    % nExperiments x 1
    % field experiments is requisite


    % make sure path has trailing filsep
    path = fullfile([path filesep]);
    if exist(path) ~= 7
        error('must supply valid path');
    end
    
    if ispref('DB', 'experiments')
        exps = getpref('DB', 'experiments');
        expIndex = find(strcmp(exps, experiment));
    else
        % initialize it
        setpref('DB', 'experiments', {});
        setpref('DB', 'paths', {});
        expIndex = [];
        disp('*** database initalized ***');
    end
    
    if expIndex        
        paths = getpref('DB', 'paths');
        paths{expIndex} = path;
        setpref('DB', 'paths', paths);
    else
        experiments = getpref('DB', 'experiments');        
        experiments{end+1} = experiment;
        setpref('DB', 'experiments', experiments);
        paths = getpref('DB', 'paths');
        paths{end+1} = path;
        setpref('DB', 'paths', paths);
    end
        
    
    
    
%     
%     % DB preference contains with fields containing cell arrays of size
%     % nExperiments x 1
%     % field experiments is requisite
%     if isempty(varargin) || rem(length(varargin), 2)
%         error('supply list of parameter value pairs via varargin');
%     end
%     
%     fieldList = {'paths'}; % right now only the path to the database file is stored
%     if ispref('DB', 'experiments')
%         exps = getpref('DB', 'experiments');
%         expIndex = find(strcmp(exps, experiment));
%     else
%         expIndex = [];
%     end
%     
%     if expIndex
%         for counter = 1:2:length(varargin
%             thisPref = getpref('Repos', 'path', path);
%             thisP
%         end
%     end