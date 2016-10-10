function peak = bpCalcPeak_dFF(Photometry, ch, window, zeroTimes, varargin)
% Photometry- output of processTrialAnalysis_Photometry, ch- channel,
% Window (size = [1,2]): start and stop time for peak calculation in seconds
% relative to zero time

% I could just have the window calculated
% This is a bare bones initial version, in future either:
% 1) have a flexible wrrapper function OR
% 2) do some of the following within this function:
% implement flexible windows (window calculation function handle?), state name(s) instead of a
% window, triggering event?
% different "peak" methods, e.g. average, max, smoothed max, fitted peak,
% etc.,  integral?

    if nargin < 4
        zeroTimes = [];
    end
    defaults = {...
        'method', 'mean';...
        'phField', 'dFF';...
%         'zeroTimes', [];... % why didn't this work as optional? isssue
%         with parseargs...
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
        zeroTimes2 = zeros(1, nTrials); % zeroTimes2 = matrix
    else
        if ~s.referenceFromEnd
            zeroTimes2 = cellfun(@(x) x(1), zeroTimes); % matrix
        else
            zeroTimes2 = cellfun(@(x) x(end), zeroTimes); % matrix
        end
    end
    for trial = 1:nTrials
        if isempty(zeroTimes) % just use window directly
            w2 = s.window;
        else
            trialStart = Photometry.startTime(trial);
            w2 = s.window + zeroTimes2(trial) - trialStart;
        end
        p1 = bpX2pnt(w2(1), Photometry.sampleRate);
        p2 = bpX2pnt(w2(2), Photometry.sampleRate);
        trialData = Photometry.data(ch).(s.phField)(trial, p1:p2);
        switch s.method
            case 'mean'
                peak.data(trial) = mean(trialData);
            case 'max'
                peak.data(trial) = max(trialData);
            otherwise
        end
    end