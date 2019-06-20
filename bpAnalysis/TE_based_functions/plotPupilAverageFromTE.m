function varargout = plotPupilAverageFromTE(TE, trials, varargin)
% vargout- first output is axis handle, second is vector of solid line handles to bounded plots

% trials- may be cell array of trials for different conditions
    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'PupilField', 'pupil';...
        'measurementField', 'pupDiameterNorm';...
        'linespec', [];...
        'window', [];...
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
    pupil = s.PupilField;
    
    if ~iscell(trials)
        trials = {trials};
    end
    
    if ~isfield(TE, pupil)
        error([pupil ' field does not exist']);
    end
    

    frameRate = max(TE.(pupil).frameRate);

    xData = TE.(pupil).xData;
    if ~isempty(s.window)
        startP = max(1, bpX2pnt(s.window(1), frameRate, xData(1)));
        endP = min(length(xData), bpX2pnt(s.window(2), frameRate, xData(1)));
    else
        startP = 1;
        endP = length(xData);
    end
    xData = xData(startP:endP);
    

    ax = s.ax;
    lh = [];
    for counter = 1:length(trials)
        thisLinespec = s.linespec{rem(counter - 1, length(s.linespec)) + 1}; % cycle through linespec if it isn't long enough        
        currentTrials = trials{counter};
        if (length(currentTrials) <= 1) || (sum(currentTrials) <= 1) % verify more than 1 trial for linear or logical indices            
            continue
        end
        currentData = TE.(pupil).(s.measurementField)(currentTrials, startP:endP);
        avg = nanmean(currentData, 1);
        avgSEM = std(currentData, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(currentData), 1));
        thisHl = boundedline(xData, avg, avgSEM, thisLinespec, ax);       
        lh(end+1) = thisHl; % return handles of the solid lines in the bounded plots
    end
    
    
    if nargout >=1
        varargout{1} = ax;
    end
    
    if nargout == 2
        varargout{2} = lh;
    end
    
    