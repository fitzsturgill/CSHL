function peak = bpCalcPeak_Pupil(pupil, window, zeroTimes, varargin)
% Pupil- output of addPupilometryToTE
% Window (size = [1,2]): start and stop time for peak calculation in seconds
% relative to zero time recording during trial

% see also bpCalcPeak_dFF

    if nargin < 4
        zeroTimes = [];
    end
    defaults = {...
        'method', 'mean';...
        'pupilField', 'dFF';...
%         'zeroTimes', [];... % why didn't this work as optional? isssue
%         with parseargs...
        'window', window;... % not optional
        'referenceFromEnd', 0;...
        };
    
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

    nTrials = size(pupil.data(ch).dFF, 1);

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
        trialStart = Photometry.startTime(trial);
        w2 = s.window + zeroTimes2(trial) - trialStart;
        p1 = bpX2pnt(w2(1), Photometry.sampleRate);
        p2 = bpX2pnt(w2(2), Photometry.sampleRate);
        trialData = Photometry.data(ch).(s.pupilField)(trial, p1:p2);
        switch s.method
            case 'mean'
                peak.data(trial) = mean(trialData);
            case 'max'
                peak.data(trial) = max(trialData);
            otherwise
        end
    end