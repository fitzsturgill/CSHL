function [ax, lh] = bpLickRaster(SessionData, type, outcome, zeroField, figName, ax)
        % create lickRaster for an individual session
        % optional arguments: if you want to use a preexisting axis,
        % then pass '' to figName and pass the axes handle
        
        % zerofield: string, e.g. 'DeliverStimulus'
    if nargin < 5
        figName = 'lickRaster';
    end
    
    if ~isempty(figName)
        fig = ensureFigure(figName, 1);
    end
    
    if nargin < 5
        clf;
    end
    
    if nargin < 6 %make a new axes unless one is provided
        ax=axes(...
        'Parent', fig,...
        'YDir', 'reverse'...
        );
    else
        axes(ax); %bring it to front and make sure that YDir is reversed
        set(ax, 'YDir', 'reverse');
    end


    [lickTimes, lickTrials, nLickTrials] = bpGetLicks(SessionData, type, outcome, zeroField);
        
    lh = linecustommarker(lickTimes, lickTrials, [], [], ax);
    % make sure that yaxis spans total number of trials so that in the absence of licks the number of lickless trials is indicated    
    try
        set(ax, 'YLim', [0, max(nLickTrials, 1)]);
    catch
        disp('wtf');
    end
        

    
    
    
    
    