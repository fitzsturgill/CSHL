function avgData = eventAverageFromTE(TE, trials, event, varargin)
% output - avgData with fields Avg, STD, SEM and N
% trials- may be cell array of trials for different conditions

% remember to include varargin-supplied parameters to
% extractEventTimesFromTE, e.g. 'zeroField' || 'zeroTimes', 'window' || 'startField' & 'endField'
% 7/2018, updated to optionally utilize zeroTime. Also updated to utilize
% window exclusive to (rather than in addition to) startField and endField
% option to use startField, endField and zeroField maintained for backward
% compatibility

    defaults = {...
        'window', [];...
        'binWidth', 0.25;... % 0.5s bins by default
        };    
    [s, ~] = parse_args(defaults, varargin{:});
    
    if ~iscell(trials)
        trials = {trials};
    end
    
    maxWindow = [min(s.window(:,1)) max(s.window(:,2))];
    
    % if interval isn't divisible by width then adjust number of bins
    binEdges = linspace(maxWindow(1), maxWindow(2), round((maxWindow(2) - maxWindow(1))/s.binWidth));
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
        Avg(counter,:) = nanmean(eventRates, 1);
        STD(counter, :) = std(eventRates, 0, 1, 'omitnan');
        SEM(counter, :) = std(eventRates,  0, 1, 'omitnan') ./ sqrt(sum(~isnan(eventRates), 1));
        N(counter, :) = sum(~isnan(eventRates), 1);
    end
    
    avgData = struct(...
        'Avg', Avg,...
        'STD', STD,...
        'SEM', SEM,...
        'N', N,...
        'xData', xData...
        );
