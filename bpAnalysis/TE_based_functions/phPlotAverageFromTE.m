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
    
    %% note I've copied zerotimes code from bpCalcPeak_dFF, this feature is not fully implemented (i.e. it just uses first element of zeroTimes
    % and currently assumes that zero time is invariant across trials
    if ~isempty(s.zeroTimes) && iscell(s.zeroTimes)
%         if ~s.referenceFromEnd
            zeroTimes2 = cellfun(@(x) x(1), s.zeroTimes); % matrix
%         else
%             zeroTimes2 = cellfun(@(x) x(end), zeroTimes); % matrix
%         end
    else
        zeroTimes2 = s.zeroTimes; % not yet tested
    end
    
%     if size(s.window, 1) == 1
%         s.window = repmat(s.window, nTrials, 1);
%     end
%%
    xData = TE.(Photometry).xData;
    if ~isempty(s.window) && ~isempty(s.zeroTimes)
        allTrialsZero = zeroTimes2(1) - TE.(Photometry).startTime(1); % Kludge adapted from bpCalcPeak_dFF
        startP = bpX2pnt(s.window(1) + allTrialsZero, TE.(Photometry).sampleRate);
        endP = bpX2pnt(s.window(2) + allTrialsZero, TE.(Photometry).sampleRate);
    elseif ~isempty(s.window) % use Photometry xData to infer zero point (really using first x value in Photometry xData)
        startP = max(1, bpX2pnt(s.window(1), TE.(Photometry).sampleRate, xData(1)));
        endP = min(length(xData), bpX2pnt(s.window(2), TE.(Photometry).sampleRate, xData(1)));        
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
    
    