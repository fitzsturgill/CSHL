function importRepositories(repos, varargin)
% recursively add repositories to matlab path, repositories correspond to
% folder names within the repository path, optional parameters-  'path',
% specify repository path
    % defaults
    if ispref('Repos', 'path')
        path = getpref('Repos', 'path');
    else
        [~, path] = uiputfile('path', 'Select default folder containing your repositories...');
        setpref('Repos', 'path', path);
    end
        
    % parse input parameter pairs
    counter = 1;
    while counter+1 <= length(varargin) 
        prop = varargin{counter};
        val = varargin{counter+1};
        switch prop
            case 'path'
                path = val;            
            otherwise
        end
        counter=counter+2;
    end
    
    if ischar(repos)
        repos = {repos};
    end
    
    for counter = 1:length(repos)
        repo = repos{counter};
        if repo(end) ~= filesep
            repo = [repo filesep];
        end
        repo = fullfile(path, repo);
        addpath(genpath(repo));
        disp(['*** added repository: ' repo ' ***']);
    end    