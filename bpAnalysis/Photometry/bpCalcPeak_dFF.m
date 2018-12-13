function peak = bpCalcPeak_dFF(Photometry, ch, window, zeroTimes, varargin)
% Photometry- output of processTrialAnalysis_Photometry, ch- channel,
% Window (size = [1,2]): start and stop time for peak calculation in seconds
% relative to zero time

% zeroTimes- can be cell array of state times (see also referenceFromEnd)
% or vector (although this option hasn't yet been tested, see
% bpCalcPeak_Pupil), must be in "photometry time" - time from photometry
% acquisition start
% zeroTimes- if omitted or empty zero is defined as beginning of photometry
% recording

    if nargin < 4
        zeroTimes = [];
    end
    % kludge to emulate zeroTimes as behaving like an optional parameter
    % (see defaults below where I chose to leave it out of varargin)
    if ischar(zeroTimes)
        varargin = horzcat({zeroTimes}, varargin);
        zeroTimes = [];
    end

    defaults = {...
        'method', 'mean';... % 'mean' or 'max', 'min', or 'percentile'
        'percentile', 0.9;... % only used with percentile calculation, default = 0.9, for 90th%
        'phField', 'dFF';...
%         'zeroTimes', [];... % why didn't this work as optional? isssue
%         with parseargs...   % note 9/2017-  I think I must have been
%         mistaken, it could work supplied as optional parameter
        'window', window;... % not optional
        'referenceFromEnd', 0;...
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

    nTrials = size(Photometry.data(ch).dFF, 1);

    peak = struct(...
        'data', zeros(nTrials, 1),...
        'settings', s,...
        'channel', ch... 
    );


    if isempty(zeroTimes)
        zeroTimes2 = Photometry.startTime; % just use start of photometry recording if zeroTimes are undefined
    elseif iscell(zeroTimes)
        if ~s.referenceFromEnd
            zeroTimes2 = cellfun(@(x) x(1), zeroTimes); % matrix
        else
            zeroTimes2 = cellfun(@(x) x(end), zeroTimes); % matrix
        end
    else
        zeroTimes2 = zeroTimes; % not yet tested
    end
    zeroTimes2 = zeroTimes2(:); 
    if isscalar(zeroTimes2)
        zeroTimes2 = repmat(zeroTimes2, nTrials, 1);
    end
        
    if size(s.window, 1) == 1
        s.window = repmat(s.window, nTrials, 1);
    end
    
    for trial = 1:nTrials
        trialZero = zeroTimes2(trial) - Photometry.startTime(trial);        
        p1 = bpX2pnt(s.window(trial,1) + trialZero, Photometry.sampleRate);
        p2 = bpX2pnt(s.window(trial,2) + trialZero, Photometry.sampleRate);
        trialData = Photometry.data(ch).(s.phField)(trial, p1:p2);
        switch s.method
            case 'mean'
                peak.data(trial) = mean(trialData);
            case 'max'
                peak.data(trial) = max(trialData);
            case 'min'
                peak.data(trial) = min(trialData);
            case 'percentile'
                peak.data(trial) = percentile(trialData, s.percentile);
            otherwise
                error('*** incorrect peak determination method ***');
        end
    end
    
    
% I could just have the window calculated
% This is a bare bones initial version, in future either:
% 1) have a flexible wrrapper function OR
% 2) do some of the following within this function:
% implement flexible windows (window calculation function handle?), state name(s) instead of a
% window, triggering event?
% different "peak" methods, e.g. average, max, smoothed max, fitted peak,
% etc.,  integral?    