function initCSHL_Fitz(rootPath)


    if nargin < 1
        [~, rootPath] = uiputfile('path', 'Choose CHL directory');
        if rootPath == 0
            return
        end
    end
    
    addpath(genpath(rootPath));
    disp('*** CSHL Fitz initialized ***');
    