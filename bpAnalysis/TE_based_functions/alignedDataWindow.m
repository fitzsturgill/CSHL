function [data, xData] = alignedDataWindow(inputData, trials, varargin)
% Fitz Sturgill 2018,  "generic" version of phAlignedWindow
% see phAlignedWindow

% returns a data array aligned to zero of size nTotalTrials x
% maxWindowSamples. 

    defaults = {...
        'zeroTimes', [];... % scalar or length TE-nTrials cell array or vector containing zero times relative to Bpod trial start        
        'window', [];... % averaging window relative to zero time/ alignment point, e.g. [-3 2] or [cellfun(@(x) x(1), TE.stateBefore) cellfun(@(x) x(end), TE.stateAfter)
        'referenceFromEnd', 0;... % relevent ONLY when zeroTimes are supplied as a cell array (e.g. if you want to align to the beginning or end of a bpod state)
        'Fs', [];...
        'startTimes', [];... % data start time relative to Bpod trial start
        };
    [s, ~] = parse_args(defaults, varargin{:});

%% argument checking    
    assert(~isempty(s.window), 'window is required parameter');
    assert(~isempty(s.zeroTimes), 'zeroTimes is required parameter');
    assert(~isempty(s.startTimes), 'data start times (in Bpod time) is a required parameter');
    assert(~isempty(s.Fs), 'Fs (sample rate) is a required parameter');
%% local variables    
    if iscell(inputData)
        uniformOutput = 0;
    else
        uniformOutput = 1;
    end
    Fs = s.Fs;

    if islogical(trials)
        trials = find(trials);
    end
    nTrials = length(trials);
    totalTrials = size(inputData, 1);
    
%% process zeroTimes and windows
    if iscell(s.zeroTimes)
        if ~s.referenceFromEnd % first point in each cell array element
            zeroTimes = cellfun(@(x) x(1), s.zeroTimes); % returns vector
        else % last point in each cell array element
            zeroTimes = cellfun(@(x) x(end), s.zeroTimes); % returns vector
        end        
    else
        zeroTimes = s.zeroTimes; 
    end    
    
    zeroTimes = zeroTimes(:); 
    assert(length(zeroTimes) == totalTrials, 'zeroTimes for all trials (typically length TE), not just the subset of trials, must be supplied');
    
    if size(s.window, 1) == 1
        s.window = repmat(s.window, totalTrials, 1);
    end
    
    if iscell(s.startTimes)
        startTimes = cellfun(@(x) x(1), s.startTimes); % returns vector
    else
        startTimes = s.startTimes;
    end
    zeroTimes = zeroTimes - startTimes; % redefine zero times relative to photometry start
    
%% initalize data array padded to maximum window size
    samplesPerWindow = ceil((max(s.window(:,2)) - min(s.window(:,1))) * Fs); 
    % NOTE! zeroPoint can be negative or > nSamples if window is DOESN'T contain zero.
    % for example, you might want a window that extends from 1-3
    % seconds after or before a zero point.
    zeroPoint = localX2pnt(0, Fs, min(s.window(:,1))); 
    xData = (0:(samplesPerWindow - 1))/Fs + min(s.window(:,1));
    
    data = NaN(nTrials, samplesPerWindow); % intialize
    
%% loop over trials to fill in data array
    for counter = 1:nTrials
        trial = trials(counter);
        zeroTime = zeroTimes(trial);
        if uniformOutput
            trialData = inputData(trial, :);
        else
            trialData = inputData{trial};
        end
        nSamples = length(trialData);
        acqDuration = nSamples / Fs;
        
        if isnan(zeroTime) || (zeroTime < 0) || (zeroTime > acqDuration)|| any(isnan(s.window(trial, :)))
            continue % continue if zero time occurs outside photometry acquisition or if windows or zeros contain NaNs
        end
        
        % source data start, end, and zero points
        szero = localX2pnt(zeroTime, Fs);
        s1 = max(localX2pnt(zeroTime + s.window(trial, 1), Fs), 1);
        s2 = min(localX2pnt(zeroTime + s.window(trial, 2), Fs), nSamples) - 1; 
        % destination data start and end points
        d1 = zeroPoint - (szero - s1);
        d2 = zeroPoint + (s2 - szero);
        % add trial data to padded array
        data(counter, d1:d2) = trialData(s1:s2);
    end
    
    
function p=localX2pnt(x, Fs, startX)

    if nargin < 3
        startX = 0;
    end
    
    deltaX=1/Fs;

    p=round(1+(x-startX)/deltaX); % local version can return negative points 
%     p=max(round(1+(x-startX)/deltaX), 1);    
        