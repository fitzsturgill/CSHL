function count = countEventFromTE(TE)
    


% remember to include varargin-supplied parameters to
% extractEventTimesFromTE, e.g. 'zeroField', 'startField', 'endField'  
    varargout = {};
    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'linespec', [];...
        'window', [];...
        'binWidth', 0.25;... % 0.5s bins by default
        'trialNumbering', 'global';... % 'singleSession'- by session, 'consecutive', 'global' cross session        };    
        };    
    [s, ~] = parse_args(defaults, varargin{:});

    if isempty(s.linespec)
        s.linespec = {'k', 'r', 'b', 'g'};
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