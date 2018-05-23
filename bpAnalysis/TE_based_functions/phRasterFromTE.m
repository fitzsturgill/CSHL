function ih = phRasterFromTE(TE, trials, ch, varargin)

    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'PhotometryField', 'Photometry';...
        'window', [];...
        'CLimFactor', 3;...  % scales color table -/+ n standard deviations away from the mean
        'trialNumbering', 'consecutive';... % 'consecutive' or 'global'
        'CLim', [];... % if specified, CLimMode set manually
        'medFilter', 0;... % if specified, set window for median filtering, only works with consecutive mode right now....
        'sortValues', [];... % overrides trialNumbering to 'consecutive', and 'showSessionBreaks', to 0
        'showSessionBreaks', 1;...
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
    
    if ~isempty(s.sortValues)
        s.showSessionBreaks = 0;
        s.trialNumbering = 'consecutive';
    end
    
    Photometry = s.PhotometryField;
    if ~isscalar(TE)
        error('TE must be scalar');
    end
    
    if ~isfield(TE, Photometry)
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
        imavg = mean(mean(TE.(Photometry).data(ch).dFF(:, startP:endP), 'omitnan'), 'omitnan');
        imstd = mean(std(TE.(Photometry).data(ch).dFF(:, startP:endP), 'omitnan'), 'omitnan');
        s.CLim = [imavg - s.CLimFactor * imstd, imavg + s.CLimFactor * imstd];
    end
    switch s.trialNumbering
        case 'consecutive'
            if ~isempty(s.sortValues)
                if islogical(trials)
                    trials = find(trials);
                end
                [sorted, key] = sort(s.sortValues(trials));
                trials = trials(key);
            end
            cData = TE.(Photometry).data(ch).dFF(trials, startP:endP);
            if s.medFilter
                cData = MEDFILT(cData, s.medFilter);
            end
            sessionBreaks = find(diff(TE.sessionIndex(trials)))';            
        %     sessionBreaks = find(diff(TE.epoch(trials)))';     % kludge for sfn poster, show epoch change (reversal)
            ih = image('Xdata', s.window, 'YData', [1 size(cData, 1)],...
                'CData', cData, 'CDataMapping', 'Scaled', 'Parent', gca);
            if s.showSessionBreaks
                line(repmat(s.window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
            end

        case 'global'
            cData = NaN(length(TE.filename), endP-startP+1);
            cData(trials,:) = TE.(Photometry).data(ch).dFF(trials, startP:endP);
            ih = image('Xdata', s.window, 'YData', [1 size(cData, 1)],...
                'CData', cData, 'CDataMapping', 'Scaled', 'Parent', gca);
            sessionBreaks = find(diff(TE.sessionIndex))';  
            if s.showSessionBreaks            
                line(repmat(s.window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks            
            end
    end
    
    set(gca, 'YLim', [1 size(cData, 1)], 'XLim', s.window, 'CLim', s.CLim);
    
end

function out = MEDFILT(cdata, window)
    out = zeros(size(cdata));
    for counter = 1:size(cdata, 1)
        out(counter, :) = medfilt1(cdata(counter, :), window);
%         out(counter, :) = smooth(cdata(counter, :), 7);
    end
end
        