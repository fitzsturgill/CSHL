function ih = phRasterFromTE(TE, trials, ch, varargin)
    
    % FS updated 7/2018, now utilizes phAlignedWindow and gains associated
    % functionality, you must pass required parameters via varargin to
    % phAlignedWindow
    % you can specify fixed or variable zeroTimes and windows, see
    % phAlignedWindow
    defaults = {...
        'ax', gca;...
        'fig', gcf;...
        'CLimFactor', 3;...  % scales color table -/+ n standard deviations away from the mean
        'trialNumbering', 'consecutive';... % 'consecutive' or 'global'
        'CLim', [];... % if specified, CLimMode set manually
        'medFilter', 0;... % if specified, set window for median filtering, only works with consecutive mode right now....
        'sortValues', [];... % overrides trialNumbering to 'consecutive', and 'showSessionBreaks', to 0
        'showSessionBreaks', 1;...
        'window', [];... % averaging window around alignment point
        'zeroTimes', [];... % empty- inferred from Photometry xData OR cell array or vector containing zero times relative to Bpod trial start        
        'PhotometryField', 'Photometry';...
        'FluorDataField', 'dFF';...
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
    
    Photometry = TE.(s.PhotometryField);
    if isempty(s.zeroTimes)
%         s.zeroTimes = valuecrossing(1:length(Photometry.xData), Photometry.xData, 0); % inferred from xData, historical usage...
%         s.zeroTimes = s.zeroTimes + Photometry.startTime; % convert to Bpod time
          s.zeroTimes = Photometry.startTime - Photometry.xData(1);
        if isempty(s.window)
            s.window = Photometry.xData([1 end]);
        end
    end
    varargin = rewrite_varargin(s, varargin);
    
    [cData, xData] = phAlignedWindow(TE, trials, ch, varargin{:});
    
    % determine CLim, use all trials and points so CLim/image scaling is consistent
    % across conditions or sets of trials (so you can compare rasters by
    % eye)
    if isempty(s.CLim)
        imavg = nanmean(nanmean(Photometry.data(ch).(s.FluorDataField), 2), 1); % FS MOD 10/2018, measure within trials and then across trials, right?
        imstd = nanmean(nanstd(Photometry.data(ch).(s.FluorDataField), 0, 2), 1);
        s.CLim = [imavg - s.CLimFactor * imstd, imavg + s.CLimFactor * imstd];
    end
    switch s.trialNumbering
        case 'consecutive'
            if ~isempty(s.sortValues)
                if islogical(trials)
                    trials = find(trials);
                end
                [~, key] = sort(s.sortValues(trials));
                trials = trials(key);
                cData = cData(key, :);
            end
            if s.medFilter
                cData = MEDFILT(cData, s.medFilter);
            end
            sessionBreaks = find(diff(TE.sessionIndex(trials)))';            
            ih = image('Xdata', [xData(1) xData(end)], 'YData', [1 size(cData, 1)],...
                'CData', cData, 'CDataMapping', 'Scaled', 'Parent', gca);
            if s.showSessionBreaks
                line(repmat(s.window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks
            end

        case 'global'
            cDataCopy = cData;
            cData = NaN(length(TE.filename), size(cDataCopy, 2));
            cData(trials,:) = cDataCopy;
            ih = image('Xdata', [xData(1) xData(end)], 'YData', [1 size(cData, 1)],...
                'CData', cData, 'CDataMapping', 'Scaled', 'Parent', gca);
            sessionBreaks = find(diff(TE.sessionIndex))';  
            if s.showSessionBreaks            
                line(repmat(s.window', 1, length(sessionBreaks)), [sessionBreaks; sessionBreaks], 'Parent', gca, 'Color', 'w', 'LineWidth', 2); % session breaks            
            end
    end    
    if size(cData, 1) > 0
        set(gca, 'YLim', [1 size(cData, 1)], 'XLim', [xData(1) xData(end)], 'CLim', s.CLim);    
    end
end

function out = MEDFILT(cdata, window)
    out = zeros(size(cdata));
    for counter = 1:size(cdata, 1)
        out(counter, :) = medfilt1(cdata(counter, :), window);
%         out(counter, :) = smooth(cdata(counter, :), 7);
    end
end
        