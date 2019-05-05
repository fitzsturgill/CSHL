function varargout = plotEventAverageFromTE(TE, trials, event, varargin)
% varargout- first output is axis handle, second is vector of solid line handles to bounded plots
% trials- may be cell array of trials for different conditions/trial subsets or vector of
% trials

% remember to include varargin-supplied parameters to
% extractEventTimesFromTE, e.g. 'zeroField' || 'zeroTimes', 'window' || 'startField' & 'endField'
% 7/2018, updated to optionally utilize zeroTime. Also updated to utilize
% window exclusive to (rather than in addition to) startField and endField
% option to use startField, endField and zeroField maintained for backward
% compatibility

varargout = {};
    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'linespec', [];...
        'window', [];...
        'binWidth', 0.25;... % 0.5s bins by default
        'cmap', [];... % for more precise control of color, supply a nTrialSets x 3 (rgb) colormap matrix
        'alpha', 1;... % 1 == transparent bounds, 0 == opaque bounds
        };    
    [s, ~] = parse_args(defaults, varargin{:});

    if isempty(s.linespec)
        s.linespec = {'g', 'r', 'b', 'k'};
    elseif ischar(s.linespec)
        s.linespec = {s.linespec};
    end
    
    if isempty(s.fig)
        s.fig = figure;
    end
    if isempty(s.ax)
        figure(s.fig);
        s.ax = axes;
    end
    
    if ~iscell(trials)
        trials = {trials};
    end
    
    if s.alpha
        alpha = {'alpha'};
    else
        alpha = {};
    end    
    

    maxWindow = [min(s.window(:,1)) max(s.window(:,2))];

    avgData = eventAverageFromTE(TE, trials, event, varargin{:});

    ax = s.ax;
    xData = avgData.xData(1,:); % fixed across trial subsets, see eventAverageFromTE
    for counter = 1:length(trials)
        thisLinespec = s.linespec{rem(counter - 1, length(s.linespec)) + 1}; % cycle through linespec if it isn't long enough     
        if any(isfinite(avgData.Avg(counter, :))) % can't be just all NaNs        
            if isempty(s.cmap)        
                thisHl = boundedline(xData, avgData.Avg(counter, :), avgData.SEM(counter, :), thisLinespec, ax, alpha{:}, 'nan', 'gap');       
            else
                thisHl = boundedline(xData, avgData.Avg(counter, :), avgData.SEM(counter, :),  ax, alpha{:}, 'cmap', s.cmap(counter, :), 'nan', 'gap');       
            end
        else
            thisHl = NaN;
        end
        lh(counter) = thisHl; % return handles of the solid lines in the bounded plots
    end
    if nargout >=1
        varargout{1} = ax;
    end
    
    if nargout == 2
        varargout{2} = lh;
    end
            
        
        
        
        
    
    

