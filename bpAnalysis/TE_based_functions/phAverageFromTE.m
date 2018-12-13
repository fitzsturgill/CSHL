function avgData = phAverageFromTE(TE, trials, ch, varargin)
% see also phAlignedWindow...
% output - avgData with fields phMean, phSTD, phSEM, N, and xData
% !!!!!! trials- may be cell array of trial subsets or a logical or linear
% index of a single trial subset
    defaults = {...
        'PhotometryField', 'Photometry';...
        'FluorDataField', 'dFF';...
        'window', [];... % averaging window around alignment point
        'zeroTimes', [];... % empty- inferred from Photometry xData OR cell array or vector containing zero times relative to Bpod trial start
        'referenceFromEnd', 0;... % relevent when zeroTimes are supplied as a cell array (e.g. if you want to align to the beginning or end of a bpod state)        
        };    
    [s, ~] = parse_args(defaults, varargin{:});

    if ~iscell(trials)
        trials = {trials};
    end
%% handle empty windows or zeroTimes and pass to phAligned window by rewriting varargin (way I'm re-writing varargin is not ideal)

    assert(isfield(TE, s.PhotometryField), [s.PhotometryField ' field does not exist']);    

    Photometry = TE.(s.PhotometryField);
    if isempty(s.zeroTimes)
        s.zeroTimes = Photometry.startTime - Photometry.xData(1);
        if isempty(s.window)
            s.window = Photometry.xData([1 end]);
        end
    end
    
    
    % rewrite varargin to pass to phAlignedWindow
    wpos = cellfun(@(x) ischar(x) && strcmp(x, 'window'), varargin);
    if ~any(wpos)
        varargin(end+1:end+2) = {'window', s.window};
    else
        varargin{find(wpos)+ 1} = s.window;
    end
    
    zpos = cellfun(@(x) ischar(x) && strcmp(x, 'zeroTimes'), varargin);
    if ~any(zpos)
        varargin(end+1:end+2) = {'zeroTimes', s.zeroTimes};
    else
        varargin{find(zpos) + 1} = s.zeroTimes;
    end
    
    

%% calculate averages
    for counter = 1:length(trials) 
        currentTrials = trials{counter};
        [data, xData] = phAlignedWindow(TE, currentTrials, ch, varargin{:});
        if counter == 1
            Avg = NaN(length(trials), length(xData));
            STD = NaN(length(trials), length(xData));
            SEM = NaN(length(trials), length(xData));
            N = NaN(length(trials), length(xData));    
            XData = NaN(length(trials), length(xData));
        end
        Avg(counter,:) = nanmean(data, 1);
        STD(counter,:) = std(data, 0, 1, 'omitnan');
        SEM(counter,:) = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));
        N(counter, :) = sum(~isnan(data), 1);
        XData(counter, :) = xData; % xData is always the same so why do I duplicate it?
    end
    
    avgData = struct(...
        'Avg', Avg,...
        'STD', STD,...
        'SEM', SEM,...
        'N', N,...
        'xData', XData... % xData is always the same so why do I duplicate it? (can't go back now...)
        );

    