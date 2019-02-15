function h = ensureFigure(figName, erase)
    % from Alex...
    % by default, do erase existing figure
    if nargin < 2
        erase = 1;
    end
    
    figs = findobj('Type', 'figure');
    names = get(figs, 'Name');
    if ischar(names)
        names = {names};
    end
    whichFigs = find(ismember(names, figName));
    
    if isempty(whichFigs)
        h = figure('Name', figName); % create figure
    else
        h = figure(figs(whichFigs(1))); %or bring to front
        if erase
            clf; % clear it for erase mode
        end
            
        if length(whichFigs > 1) % and close any duplicates
            close(figs(whichFigs(2:end)));
        end
    end



