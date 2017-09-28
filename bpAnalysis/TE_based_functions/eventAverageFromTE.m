function avgData = eventAverageFromTE(TE, trials, event, varargin)
% output - avgData with fields Avg, STD, SEM, and N
% trials- may be cell array of trials for different conditions

% remember to include varargin-supplied parameters to
% extractEventTimesFromTE, e.g. 'zeroField', 'startField', 'endField'  
    defaults = {...
        'window', [];...
        'binWidth', 0.25;... % 0.5s bins by default
        };    
    [s, ~] = parse_args(defaults, varargin{:});
    
    if ~iscell(trials)
        trials = {trials};
    end    
    
    binEdges = linspace(s.window(1), s.window(2), round((s.window(2) - s.window(1))/s.binWidth)); % if interval isn't divisible by width then adjust number of bins
    s.binWidth = binEdges(2) - binEdges(1); % and adjust width if necessary
    % initialize 
    Avg = NaN(length(trials), length(binEdges) - 1);
    STD = NaN(length(trials), length(binEdges) - 1);
    SEM = NaN(length(trials), length(binEdges) - 1);    
    N = NaN(length(trials), length(binEdges) - 1);
    
    binCenters = binEdges(1:end-1) + s.binWidth;
    xData = repmat(binCenters, length(trials), 1);

    for counter = 1:length(trials)
        currentTrials = trials{counter};
        nTrials = length(find(currentTrials));   % currentTrials can be logical or linear indexes     
        [eventTimes, eventTrials] = extractEventTimesFromTE(TE, currentTrials, event, varargin{:}); 
        counts = histCountsByTrial(eventTimes, eventTrials, binEdges);
        eventRates = counts / s.binWidth;
        Avg(counter,:) = mean(eventRates);
        STD(counter, :) = std(eventRates);
        SEM(counter, :) = std(eventRates) ./ sqrt(nTrials);
        N(counter, :) = nTrials;
    end
    
    avgData = struct(...
        'Avg', Avg,...
        'STD', STD,...
        'SEM', SEM,...
        'N', N,...
        'xData', xData...
        );
