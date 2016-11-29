function ih = phRasterFromTE(TE, trials, ch, varargin)

    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'PhotometryField', 'Photometry';...
        'window', [];...
        'CLimFactor', 3;...
        'CLim', [];... % if specified, CLimMode set manually
        };
    [s, ~] = parse_args(defaults, varargin{:});
    if isempty(s.fig)
        s.fig = figure;
    else
        figure(s.fig);
    end
    
    if isempty(s.ax)
        figure(s.fig);
        s.ax = axes('YDir', 'Reverse');
    else
        set(s.ax, 'YDir', 'Reverse');
    end
    
    Photometry = s.PhotometryField;
    if ~isscalar(TE)
        error('TE must be scalar');
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
        s.window = [xData(1) xData(end)];
    end
    
    

    % determine CLim, use all trials so CLim/image scaling is consistent
    % across conditions or sets of trials
    if isempty(s.CLim)
        imavg = mean(mean(TE.(Photometry).data(ch).dFF(:, startP:endP), 'omitnan'));
        imstd = mean(std(TE.(Photometry).data(ch).dFF(:, startP:endP), 'omitnan'));
        s.CLim = [imavg - s.CLimFactor * imstd, imavg + s.CLimFactor * imstd];
    end
    cData = TE.(Photometry).data(ch).dFF(trials, startP:endP);
    
    
    sessionBreaks = find(diff(TE.sessionIndex(trials)))';            
%     sessionBreaks = find(diff(TE.epoch(trials)))';     % kludge for sfn poster, show epoch change (reversal)
    ih = image('Xdata', s.window, 'YData', [1 size(cData, 1)],...
        'CData', cData, 'CDataMapping', 'Scaled', 'Parent', gca);
    line(repmat(s.window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
    set(gca, 'YLim', [1 size(cData, 1)], 'XLim', s.window, 'CLim', s.CLim);