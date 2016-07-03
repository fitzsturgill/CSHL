function sessions = bpLoadSession(filename, filepath)
    % loads a single session
    % creates a folder at level of that session for analysis
    % cd into that folder
    
    %arguments
    %filename- string
    %filepath - string
    
    if nargin < 2
        filepath = pwd;
    end
    
    cd(filepath);
    if nargin < 1 || isempty(filename) % interactive loading
        [filename, filepath] = uigetfile('*Session*.mat', 'select session');
        if ~filename
            sessions=[];
            return
        end
    end
    
    sessions = bpLoadSessions([], filename, filepath); % load the session
    [~, filename, ~] = fileparts(filename);
    newFolder = [filepath filename];
    if ~isdir(newFolder);
        mkdir([newFolder '\']);
    end
    cd(newFolder);
    