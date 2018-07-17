function varargout = phPlotAverageFromTE(TE, trials, ch, varargin)
% vargout- first output is axis handle, second is vector of solid line handles to bounded plots

% trials- may be cell array of trials for different conditions
    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'PhotometryField', 'Photometry';...
        'FluorDataField', 'dFF';...
        'linespec', [];... % {'k', 'r'}
        'window', [];... % window to plot with respect to zero time that is already calculated by processTrialAnalysis_Photometry2 (I think !!!!)
        'zeroTimes', [];... % not fully implemented, see usage below
        'cmap', [];...
        'alpha', 1;...  % you can
        'referenceFromEnd', 0;... % relevent when zeroTimes are supplied as a cell array (e.g. if you want to align to the beginning or end of a bpod state)
        };    
    [s, ~] = parse_args(defaults, varargin{:});


    if isempty(s.linespec)
        s.linespec = {'k', 'r', 'b', 'g'};
    elseif ischar(s.linespec) && strcmp(s.linespec, 'default')
        s.linespec = {''}; % use figure/axes colormap
    elseif ischar(s.linespec)
        s.linespec = {s.linespec};
    end
    if isempty(s.fig)
        s.fig = figure;
    end
    if isempty(s.ax)
        figure(s.fig);
        s.ax = axes;
    end
    
    if s.alpha
        alpha = {'alpha'};
    else
        alpha = {};
    end
    
    
    if ~iscell(trials)
        trials = {trials};
    end
    
    if ~isfield(TE, s.PhotometryField)
        error([Photometry ' field does not exist']);
    else
        Photometry = TE.(s.PhotometryField);
    end
    

    xData = Photometry.xData;
    if ~isempty(s.window) && ~isempty(s.zeroTimes)
        allTrialsZero = zeroTimes2(1) - Photometry.startTime(1); % Kludge adapted from bpCalcPeak_dFF
        startP = bpX2pnt(s.window(1) + allTrialsZero, Photometry.sampleRate);
        endP = bpX2pnt(s.window(2) + allTrialsZero, Photometry.sampleRate);
    elseif ~isempty(s.window) % use Photometry xData to infer zero point (really using first x value in Photometry xData)
        startP = max(1, bpX2pnt(s.window(1), Photometry.sampleRate, xData(1)));
        endP = min(length(xData), bpX2pnt(s.window(2), Photometry.sampleRate, xData(1)));        
    else
        startP = 1;
        endP = length(xData);
        s.window = [xData(1) xData(end)];
    end
%     xData = xData(startP:endP);
%% rewrite xData, Kludge, because it doesn't fully implement zeroTimes as above
    xData = linspace(s.window(1), s.window(2), endP - startP + 1);  % rewrite xData
%%
    

    ax = s.ax;
    for counter = 1:length(trials)
        thisLinespec = s.linespec{rem(counter - 1, length(s.linespec)) + 1}; % cycle through linespec if it isn't long enough        
        currentTrials = trials{counter};
        currentData = Photometry.data(ch).(s.FluorDataField)(currentTrials, startP:endP);
        avg = nanmean(currentData);
        avgSEM = std(currentData, 'omitnan') ./ sqrt(sum(~isnan(currentData), 1));

        if isempty(s.cmap)
            [thisHl, thisHp] = boundedline(xData, avg, avgSEM, thisLinespec, ax, alpha{:});       
        else
            [thisHl, thisHp] = boundedline(xData, avg, avgSEM, ax, alpha{:}, 'cmap', s.cmap(counter,:));    
        end
        lh(counter) = thisHl; % return handles of the solid lines in the bounded plots
        ph(counter) = thisHp;
    end
    
    
    if nargout >=1
        varargout{1} = ax;
    end
    
    if nargout >= 2
        varargout{2} = lh;
    end
    
    if nargout == 3
        varargout{3} = ph;
    end
    
    