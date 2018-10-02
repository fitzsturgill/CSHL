function success = dbLoadAnimal(DB, animal, varargin)
% loads experiment variables into calling workspace 
% db-   experiment database
    
    success = 1;
    try
        evalin('caller', sprintf('load(''%s'')', fullfile(DB.path, 'animals', animal, 'TE.mat')));
        if isfield(DB, 'conditions') && ~isempty(DB.conditions)            
            evalin('caller', DB.conditions);
        end
    catch ME
        ME.message
        success = 0;
    end 