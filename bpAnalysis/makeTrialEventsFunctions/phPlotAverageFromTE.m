function varargout = phPlotAverageFromTE(TE, trials, ch, varargin)
% vargout- first output is axis handle, second is vector of solid line handles to bounded plots

% trials- may be cell array of trials for different conditions
    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'PhotometryField', 'Photometry';...
        'FluorDataField', 'dFF';...
        'linespec', [];...
        'window', [];...
        };    
    [s, ~] = parse_args(defaults, varargin{:});

%     if isempty(s.PhotometryField)
%         s.PhotometryField = 'Photometry';
%     end
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
    Photometry = s.PhotometryField;
    
    if ~iscell(trials)
        trials = {trials};
    end
    
    if ~isfield(TE, Photometry);
        error([Photometry ' field does not exist']);
    end
    
    xData = TE.(Photometry).xData;
    if ~isempty(s.window)
        startP = max(1, bpX2pnt(s.window(1), TE.(Photometry).sampleRate, xData(1)));
        endP = min(length(xData), bpX2pnt(s.window(2), TE.(Photometry).sampleRate, xData(1)));
    else
        startP = 1;
        endP = length(xData);
    end
    xData = xData(startP:endP);
    

    ax = s.ax;
    for counter = 1:length(trials)
        thisLinespec = s.linespec{rem(counter - 1, length(s.linespec)) + 1}; % cycle through linespec if it isn't long enough        
        currentTrials = trials{counter};
        currentData = TE.(Photometry).data(ch).(s.FluorDataField)(currentTrials, startP:endP);
        avg = nanmean(currentData);
        avgSEM = std(currentData, 'omitnan') ./ sqrt(sum(~isnan(currentData), 1));
        thisHl = boundedline(xData, avg, avgSEM, thisLinespec, ax, 'alpha');       
        lh(counter) = thisHl; % return handles of the solid lines in the bounded plots
    end
    
    
    if nargout >=1
        varargout{1} = ax;
    end
    
    if nargout == 2
        varargout{2} = lh;
    end
    
    