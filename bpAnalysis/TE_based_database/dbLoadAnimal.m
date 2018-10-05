function success = dbLoadAnimal(DB, animal, varargin)
% loads experiment variables into calling workspace 
% db-   experiment database

    %% optional parameters, first set defaults
    defaults = {...
        'loadConditions', 1;...            % load conditions script (creates trial subset variables in caller worksapce)?
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings     
    
    success = 1;
    try
        evalin('caller', sprintf('load(''%s'')', fullfile(DB.path, 'animals', animal, 'TE.mat')));
        if s.loadConditions && isfield(DB, 'conditions') && ~isempty(DB.conditions)            
            evalin('caller', DB.conditions);
        end
    catch ME
        ME.message
        success = 0;
    end