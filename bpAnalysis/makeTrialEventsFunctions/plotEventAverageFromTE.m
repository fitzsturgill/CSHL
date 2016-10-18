function varargout = plotEventAverageFromTE(TE, trials, event, varargin)
% varargout- first output is axis handle, second is vector of solid line handles to bounded plots
% trials- may be cell array of trials for different conditions

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
    

    
    binEdges = linspace(s.window(1), s.window(2), round((s.window(2) - s.window(1))/s.binWidth)); % if interval isn't divisible by width then adjust number of bins
    s.binWidth = binEdges(2) - binEdges(1); % and adjust width if necessary
    % initialize 
    lickRates = zeros(length(trials), length(binEdges) - 1);
    lickRatesN = zeros(1, length(binEdges) - 1);
    binCenters = binEdges(1:end-1) + s.binWidth;
    

    ax = s.ax;
    for counter = 1:length(trials)
        thisLinespec = s.linespec{rem(counter - 1, length(s.linespec)) + 1}; % cycle through linespec if it isn't long enough        
        currentTrials = trials{counter};        
        nTrials = length(find(currentTrials));        
        [eventTimes, ~] = extractEventTimesFromTE(TE, currentTrials, event, varargin{:}); 
        [counts, ~] = histcounts(eventTimes, binEdges);
        eventRates = counts /nTrials /s.binWidth;
        eventRatesSEM = std(eventRates) ./ sqrt(nTrials);
        thisHl = boundedline(binCenters, eventRates, eventRatesSEM, thisLinespec, ax, 'alpha');       
        lh(counter) = thisHl; % return handles of the solid lines in the bounded plots
    end
    if nargout >=1
        varargout{1} = ax;
    end
    
    if nargout == 2
        varargout{2} = lh;
    end
            
        
        
        
        
    
    

