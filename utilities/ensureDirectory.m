function ensureDirectory(fullpath)
% Fitz Sturgill 2016
% recursively creates folders along a directory path
% fullpath = fully specified directory path


    subpath = '';
    while 1
        [fname, fullpath] = strtok(fullpath, filesep);
        if isempty(fname)
            break
        elseif ~exist(fname, 'dir')
            if isempty(subpath)
                error('root folder needs to exist');
            end
            mkdir(subpath, fname);
        end
        subpath = fullfile(subpath, fname);
    end
    
    
        