function success = dbInitExperiment(varargin)
% initializes a database

%% optional parameters, first set defaults
defaults = {...
    'name', '';...
    'conditions', '';... % name of script to generate sets of trials (to easily recycle code from my analysis scripts which operate in the base workspace)
    'animals', [];... % will contain names/folder specifiers of animal subjects
    };
[s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
success = 1;
basepath = uigetdir; % prompts windows/mac osx to give you a location to save
if isempty(basepath)
    success = 0;
    return
end

if isempty(s.name)
    % just use folder name if empty
    seps = basepath == filesep;
    sep = find(seps(1:end-1), 1, 'last');
    s.name = basepath(sep+1:end);
    if s.name(end) == filesep
        s.name = s.name(1:end-1);
    end
end

ensureDirectory(fullfile(basepath, 'animals', filesep));
ensureDirectory(fullfile(basepath, 'pooled', filesep));

DB = struct(...
    'path', basepath...
    );

save(fullfile(basepath, 'DB.mat'), 'DB');

dbRegisterExperiment(s.name, basepath);
disp(['*** Saved: ' fullfile(basepath, 'DB.mat')]);