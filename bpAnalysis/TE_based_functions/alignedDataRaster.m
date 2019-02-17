function alignedDataRaster(inputData, trials, varargin)
% Fitz Sturgill 2/2019
% see phRasterFromTE and alignedDataWindow
% inputData: matrix or cell array of height n total trials (first
% dimension)
defaults = {...
    'startTimes', [];... % data start time relative to Bpod trial start    
    'zeroTimes', [];... % zero times relative to Bpod trial start, scalar or length TE-nTrials cell array or vector containing         
    'window', [];... % averaging window relative to zero time/ alignment point, e.g. [-3 2] or [cellfun(@(x) x(1), TE.stateBefore) cellfun(@(x) x(end), TE.stateAfter)
    'referenceFromEnd', 0;... % relevent ONLY when zeroTimes are supplied as a cell array (e.g. if you want to align to the beginning or end of a bpod state)
    'Fs', [];...
    'ax', gca;...
    'fig', gcf;...
    'CLimFactor', 3;...  % scales color table -/+ n standard deviations away from the mean
    'trialNumbering', 'consecutive';... % 'consecutive' or 'global'
    'CLim', [];... % if specified, CLimMode set manually
    'sortValues', [];... % overrides trialNumbering to 'consecutive'
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
totalTrials = size(inputData, 1);



if isempty(s.CLim)
    if isnumeric(inputData)
        imavg = nanmean(nanmean(inputData, 2), 1); 
        imstd = nanmean(nanstd(inputData, 0, 2), 1);            
    elseif iscell(inputData)
        imavg = nanmean(cellfun(@(x) nanmean(x, 2), inputData), 1); 
        imstd = nanmean(cellfun(@(x) nanstd(x, 0, 2), inputData), 1);                        
    end
    s.CLim = [imavg - s.CLimFactor * imstd, imavg + s.CLimFactor * imstd];
end

if size(s.window, 1) == 1
    s.window = repmat(s.window, totalTrials, 1);
end

varargin = rewrite_varargin(s, varargin);
[cData, xData] = alignedDataWindow(inputData, trials, varargin{:});

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
        ih = image('Xdata', [xData(1) xData(end)], 'YData', [1 size(cData, 1)],...
            'CData', cData, 'CDataMapping', 'Scaled', 'Parent', gca);
    case 'global'
        cDataCopy = cData;
        cData = NaN(totalTrials, size(cDataCopy, 2));
        cData(trials,:) = cDataCopy;
        ih = image('Xdata', [xData(1) xData(end)], 'YData', [1 size(cData, 1)],...
            'CData', cData, 'CDataMapping', 'Scaled', 'Parent', gca);
    otherwise
end

if size(cData, 1) > 0
    set(gca, 'YLim', [1 size(cData, 1)], 'XLim', [xData(1) xData(end)], 'CLim', s.CLim);    
end



