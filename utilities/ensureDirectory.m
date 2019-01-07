function fullpath = ensureDirectory(fullpath)
% Fitz Sturgill 2016
% recursively creates folders along a directory path
% fullpath = fully specified directory path
    
    if ispc
        subpath = '';
    else
        subpath = '/'; % assumed to be a mac
    end
    while 1
        [fname, fullpath] = strtok(fullpath, filesep);     
        if isempty(fname)
            break
        elseif ~exist(fullfile(subpath, fname), 'dir')
            mkdir(subpath, fname);
        end
        subpath = fullfile(subpath, fname);           
    end