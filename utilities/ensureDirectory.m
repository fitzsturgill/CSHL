function ensureDirectory(fullpath)
% Fitz Sturgill 2016
% recursively creates folders along a directory path
% fullpath = fully specified directory path


    subpath = '';
    while 1
        [fname, fullpath] = strtok(fullpath, filesep);
        subpath = fullfile(subpath, fname);        
        if isempty(fname)
            break
        elseif ~exist(subpath, 'dir')
            mkdir(subpath, fname);
        end
    end
    
    
        