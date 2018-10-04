function [data, xData] = phAlignedWindow(TE, trials, ch, varargin)
% returns a photometry data array aligned to zero of size nTotalTrials x
% maxWindowSamples. Zeros for each and every trial are defined relative to Bpod trial start and may
% be supplied as vectors or cell arrays of length nTotalTrials (in which case you can select first
% or last element to which to align)
% totalTrials -> length of TE.filename not count of subset of trials passed
% via 'trials' argument
    defaults = {...
        'PhotometryField', 'Photometry';...
        'FluorDataField', 'dFF';...
        'zeroTimes', [];... % scalar or length TE-nTrials cell array or vector containing zero times relative to Bpod trial start        
        'window', [];... % averaging window relative to zero time/ alignment point, e.g. [-3 2] or [cellfun(@(x) x(1), TE.stateBefore) cellfun(@(x) x(end), TE.stateAfter)
        'referenceFromEnd', 0;... % relevent ONLY when zeroTimes are supplied as a cell array (e.g. if you want to align to the beginning or end of a bpod state)        
        };
    [s, ~] = parse_args(defaults, varargin{:});

%% argument checking    
    assert(~isempty(s.window), 'window is required parameter');
    assert(~isempty(s.zeroTimes), 'zeroTimes is required parameter');
    assert(isfield(TE, s.PhotometryField), [s.PhotometryField ' field does not exist']);
    
%% local variables    
    Photometry = TE.(s.PhotometryField);   
    Fs = Photometry.sampleRate;

    if islogical(trials)
        trials = find(trials);
    end
    nTrials = length(trials);
    totalTrials = length(TE.filename);
    
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
    assert(length(zeroTimes) == totalTrials, 'zeroTimes for entire TE (not just the subset of trials) must be supplied');
        
    if size(s.window, 1) == 1
        s.window = repmat(s.window, totalTrials, 1);
    end
    
    zeroTimes = zeroTimes - TE.Photometry.startTime; % redefine zero times relative to photometry start
    

%% initalize data array padded to maximum window size
    if Photometry.settings.uniformOutput
        samplesPerWindow = ceil((max(s.window(:,2)) - min(s.window(:,1))) * Fs); 
        % NOTE! zeroPoint can be negative or > nSamples if window is DOESN'T contain zero.
        % for example, you might want a window that extends from 1-3
        % seconds after or before a zero point.
        zeroPoint = localX2pnt(0, Fs, min(s.window(:,1))); 
        xData = (0:(samplesPerWindow - 1))/Fs + min(s.window(:,1));
    else
        % FINISH CODING!!!!
        error('finish coding to allow different sized photometry acquisitions');
    end
    
    data = NaN(nTrials, samplesPerWindow); % intialize
    
%% loop over trials to fill in data array
    for counter = 1:nTrials
        trial = trials(counter);
        zeroTime = zeroTimes(trial);
        if Photometry.settings.uniformOutput
            trialData = Photometry.data(ch).(s.FluorDataField)(trial, :);
        else
            % FINISH CODING!!!!
            % trialData = Photometry.data(ch).(s.FluorDataField){trial};???
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
        