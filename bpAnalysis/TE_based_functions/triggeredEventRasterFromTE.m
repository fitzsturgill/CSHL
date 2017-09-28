function [ax, lh] = triggeredEventRasterFromTE(TE, event, TS, varargin)
    
    % remember to include varargin-supplied parameters to
    % extractEventTimesFromTE- 'zeroField', 'startField', 'endField'
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
    tre = extractTriggeredEvents(TE, event, TS, varargin{:});
    lh = linecustommarker(tre.eventTimes, tre.eventStampIndices, [], [], s.ax);
    ax = s.ax;
    