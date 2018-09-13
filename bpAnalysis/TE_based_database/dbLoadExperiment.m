function success = dbLoadExperiment(db, animal, varargin)
% loads experiment variables into calling workspace 
% db-   experiment database
    
    success = 1;
    try
        evalin('caller', sprintf('load(''%s'')', fullfile(db.path, 'animals', animal, 'TE.mat')));
        if isfield(db, 'conditions') && ~isempty(db.conditions)            
            evalin('caller', db.conditions);
        end
    catch ME
        ME.message
        success = 0;
    end 