function success = dbSaveAnimal(DB, animal)

    % save TE
    success = 1;
%     try
        filepath = fullfile(DB.path, 'animals', animal, 'TE.mat');
        evalin('caller', sprintf('save(''%s'', ''TE'')', filepath));
        fprintf('*** saved %s ***\n', filepath);
%     catch ME
%         ME.message
%         success = 0;
%     end