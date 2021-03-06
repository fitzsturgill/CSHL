function varargout = phPlotAverageFromTE(TE, trials, ch, varargin)
% vargout- first output is axis handle, second is vector of solid line handles to bounded plots

% trials- may be cell array of trials for different conditions
    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'PhotometryField', 'Photometry';...
        'FluorDataField', 'dFF';...
        'linespec', [];... % {'k', 'r'}
        'window', [];... 
        'zeroTimes', [];... 
        'cmap', [];...
        'alpha', 1;...  % 1 == transparent bounds, 0 == opaque bounds
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
    
    assert(isfield(TE, s.PhotometryField), [s.PhotometryField ' field does not exist']);
    Photometry = TE.(s.PhotometryField);
    
    avgData = phAverageFromTE(TE, trials, ch, varargin{:}); % compile averages

    ax = s.ax;
    nLines = size(avgData.Avg, 1);
    xData = avgData.xData(1,:);    % fixed across trial subsets, see phAverageFromTE
    for counter = 1:nLines
        thisLinespec = s.linespec{rem(counter - 1, length(s.linespec)) + 1}; % cycle through linespec if it isn't long enough        
        if any(isfinite(avgData.Avg(counter, :))) % can't be just all NaNs
            if isempty(s.cmap)
                [thisHl, thisHp] = boundedline(xData, avgData.Avg(counter, :), avgData.SEM(counter, :), thisLinespec, ax, alpha{:}, 'nan', 'gap');       
            else
                [thisHl, thisHp] = boundedline(xData, avgData.Avg(counter, :), avgData.SEM(counter, :), ax, alpha{:}, 'cmap', s.cmap(counter,:), 'nan', 'gap');    
            end
            thisHl = thisHl;
            thisHp = thisHp;
        else
            continue % stupid matlab now handles line and patch handles as options, yet doesn't allow you to generate null handles with NaNs
%             thisHl = NaN;
%             thisHp = NaN;
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
    
    