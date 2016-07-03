function phPlotSessionAvg(SessionData, type, outcome, figName, ax, linespec)
    
    if nargin < 6
        linespec = {'k', 'r', 'b', 'g'};
    elseif ischar(linespec)
        linspec = {linespec};
    end

    if nargin < 4
        figName = 'photometryAvg';
    end
    
    if ~isempty(figName)
        fig = ensureFigure(figName);
    end
    
    if nargin < 4
        clf
    end
    
    if nargin < 5 || isempty(ax) %make a new axes unless one is provided
        ax=axes(...
        'Parent', fig...
        );
    else
        axes(ax);
    end

    for counter = 1:length(type)
        thisType = type(counter);
        thisOutcome = outcome(counter);
        thisLinespec = linespec{counter};
        phSessionAvg(SessionData, thisType, thisOutcome, ax, thisLinespec);
    end
        
        
        
        
        
        
        
        
        
        