function varargout = phAverageFromTE(TE, trials, ch, varargin)
% vargout- first output is axis handle, second is vector of solid line handles to bounded plots

% trials- may be cell array of trials for different conditions
    defaults = {...
        'PhotometryField', 'Photometry';...
        'window', [];...
        };    
    [s, ~] = parse_args(defaults, varargin{:});


    Photometry = s.PhotometryField;    

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
    



    currentTrials = trials; % vestigal
    currentData = TE.(Photometry).data(ch).dFF(currentTrials, startP:endP);
    avg = nanmean(currentData);
    avgSEM = std(currentData, 'omitnan') ./ sqrt(sum(~isnan(currentData), 1));

    
    
    if nargout >=1
        varargout{1} = avg;
    end
    
    if nargout == 2
        varargout{2} = avgSEM;
    end
    
    