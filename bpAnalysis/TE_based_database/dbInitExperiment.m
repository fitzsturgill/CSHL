function success = dbInitExperiment(varargin)
% initializes a database

%% optional parameters, first set defaults
defaults = {...
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

ensureDirectory(fullfile(basepath, 'animals', filesep));
ensureDirectory(fullfile(basepath, 'pooled', filesep));

DB = struct(...
    'path', basepath...
    );

save(fullfile(basepath, 'DB.mat'), 'DB');
disp(['*** Saved: ' fullfile(basepath, 'DB.mat')]);