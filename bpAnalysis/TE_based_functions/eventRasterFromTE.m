function [ax, lh] = eventRasterFromTE(TE, trials, event, varargin)
    
    % remember to include varargin-supplied parameters to
% extractEventTimesFromTE, e.g. 'zeroField' || 'zeroTimes', 'window' || 'startField' & 'endField'
    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'trialNumbering', 'global';... % 'singleSession'- by sesssion, 'consecutive', 'global' cross session        };
        'sortValues', [];... % overrides trial numbering to 'consecutive'
        'LineWidth', 0.75;...
        };
    [s, ~] = parse_args(defaults, varargin{:});
    if isempty(s.fig)
        s.fig = figure;
    end
    if isempty(s.ax)
        figure(s.fig);
        s.ax = axes;
    end
    if ~isempty(s.sortValues)
        s.trialNumbering = 'consecutive';
        if islogical(trials)
            trials = find(trials);
        end
        [~, key] = sort(s.sortValues(trials));
        trials = trials(key);
    end
    set(s.ax, 'YDir', 'reverse');

    [eventTimes, eventTrials] = extractEventTimesFromTE(TE, trials, event, varargin{:});
    lh = linecustommarker(eventTimes, eventTrials, [], [], s.ax);
    set(lh, 'LineWidth', s.LineWidth, 'Color', [0 0 0]);
    ax = s.ax;
    switch s.trialNumbering
        case 'consecutive'
            if islogical(trials) && (sum(trials) > 0)
                set(s.ax, 'YLim', [0 sum(trials)]);
            elseif ~isempty(trials)
                set(s.ax, 'YLim', [0 length(trials)]);
            end
        case 'global'
            set(s.ax, 'YLim', [0 length(TE.filename)]);
    end
    