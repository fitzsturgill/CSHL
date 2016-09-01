function [ax, lh] = eventRasterFromTE(TE, trials, event, varargin)
    
    % remember to include varargin-supplied parameters to
    % extractEventTimesFromTE
    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'trialNumbering', 'global';... % 'singleSession'- by sesssion, 'consecutive', 'global' cross session        };
        };
    [s, ~] = parse_args(defaults, varargin{:});
    if isempty(s.fig)
        s.fig = figure;
    end
    if isempty(s.ax)
        figure(s.fig);
        s.ax = axes;
    end
    set(s.ax, 'YDir', 'reverse');
    [eventTimes, eventTrials] = extractEventTimesFromTE(TE, trials, event, varargin{:});
    lh = linecustommarker(eventTimes, eventTrials, [], [], s.ax);
    ax = s.ax;
    