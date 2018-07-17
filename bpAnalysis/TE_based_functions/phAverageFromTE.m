function avgData = phAverageFromTE(TE, trials, ch, varargin)
% output - avgData with fields phMean, phSTD, phSEM, and N

% trials- may be cell array of trials for different conditions
    defaults = {...
        'PhotometryField', 'Photometry';...
        'FluorDataField', 'dFF';...
        'window', [];... % window to plot with respect to zero time that is already calculated by processTrialAnalysis_Photometry2 (I think !!!!)
        'zeroTimes', [];... % scalar, assumes consistent trial timing
        'referenceFromEnd', 0;... % relevent when zeroTimes are supplied as a cell array (e.g. if you want to align to the beginning or end of a bpod state)        
        };    
    [s, ~] = parse_args(defaults, varargin{:});


    Photometry = s.PhotometryField;
    
    if ~iscell(trials)
        trials = {trials};
    end
    
    if ~isfield(TE, Photometry)
        error([Photometry ' field does not exist']);
    end
    
%%
    if isempty(s.zeroTimes)
        zeroTimes = valuecrossing(1:length(TE.Photometry.xData), TE.Photometry.xData, 0); % inferred from xData, historical usage...
        if isempty(s.window)
            s.window = TE.Photometry.xData([1 end]);
        end
    elseif iscell(s.zeroTimes)
        if ~s.referenceFromEnd
            zeroTimes = cellfun(@(x) x(1), s.zeroTimes); % returns vector
        else
            zeroTimes = cellfun(@(x) x(end), zeroTimes); % returns vector
        end
        assert(~isempty(s.window), 'user must supply window if zero time(s) are supplied as well');
    else
        zeroTimes = s.zeroTimes; 
        assert(~isempty(s.window), 'user must supply window if zero time(s) are supplied as well');
    end
    zeroTimes = zeroTimes(:); 
    if isscalar(zeroTimes)
        zeroTimes = repmat(zeroTimes, nTrials, 1);
    end
        
    if size(s.window, 1) == 1
        s.window = repmat(s.window, nTrials, 1);
    end
    

%% determine maximum number of data points preceding and initalize 

    
    if Photometry.settings.uniformOutput
    
    end

    xData = TE.(Photometry).xData;
    if ~isempty(s.window) && ~isempty(s.zeroTimes)
        allTrialsZero = zeroTimes(1) - TE.(Photometry).startTime(1);
        startP = max(1, bpX2pnt(s.window(1) + allTrialsZero, TE.(Photometry).sampleRate));
        endP = min(length(xData), bpX2pnt(s.window(2) + allTrialsZero, TE.(Photometry).sampleRate));
    elseif ~isempty(s.window) % use Photometry xData to infer zero point (really using first x value in Photometry xData)
        startP = max(1, bpX2pnt(s.window(1), TE.(Photometry).sampleRate, xData(1)));
        endP = min(length(xData), bpX2pnt(s.window(2), TE.(Photometry).sampleRate, xData(1)));        
    else
        startP = 1;
        endP = length(xData);
        s.window = [xData(1) xData(end)];
    end

%% rewrite xData
    xData = linspace(s.window(1), s.window(2), endP - startP + 1);  % rewrite xData
    
    Avg = NaN(length(trials), length(xData));
    STD = NaN(length(trials), length(xData));
    SEM = NaN(length(trials), length(xData));
    N = NaN(length(trials), length(xData));    
    XData = NaN(length(trials), length(xData));        
    for counter = 1:length(trials)        
        currentTrials = trials{counter};
        currentData = TE.(Photometry).data(ch).(s.FluorDataField)(currentTrials, startP:endP);
        Avg(counter,:) = nanmean(currentData);
        STD(counter,:) = std(currentData, 'omitnan');
        SEM(counter,:) = std(currentData, 'omitnan') ./ sqrt(sum(~isnan(currentData), 1));
        N(counter, :) = sum(~isnan(currentData), 1);
        XData(counter, :) = xData;
    end
    
    avgData = struct(...
        'Avg', Avg,...
        'STD', STD,...
        'SEM', SEM,...
        'N', N,...
        'xData', XData...
        );

    