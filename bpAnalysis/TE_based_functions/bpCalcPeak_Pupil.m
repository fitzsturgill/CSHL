function peak = bpCalcPeak_Pupil(pupil, window, zeroTimes, varargin)
% Pupil- output of addPupilometryToTE
% Window (size = [1,2]): start and stop time for peak calculation in seconds
% relative to zero time recording during trial
% zeroTimes- can be cell array of state times (see also referenceFromEnd)
% or vector (although this option hasn't yet been tested)

% see also bpCalcPeak_dFF

    if nargin < 4
        zeroTimes = [];
    end
    defaults = {...
        'method', 'mean';... % mean or peak
        'pupilField', 'pupDiameterNorm';...
%         'zeroTimes', [];... % why didn't this work as optional? isssue
%         with parseargs...
        'window', window;... % not optional
        'referenceFromEnd', 0;...
        };
    
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

    nTrials = size(pupil.startTime, 1);

    peak = struct(...
        'data', zeros(nTrials, 1),...
        'settings', s...
    );


    if isempty(zeroTimes)
        zeroTimes2 = pupil.startTime; % use start of pupil recording if zeroTimes are undefined
    elseif iscell(zeroTimes)
        if ~s.referenceFromEnd
            zeroTimes2 = cellfun(@(x) x(1), zeroTimes); % matrix
        else
            zeroTimes2 = cellfun(@(x) x(end), zeroTimes); % matrix
        end
    else
        zeroTimes2 = zeroTimes; % not yet tested
    end
    
    for trial = 1:nTrials
        frameRate = pupil.frameRate(trial);        
        nFrames = length(pupil.(s.pupilField)(trial, :));
        trialZero = zeroTimes2(trial) - pupil.startTime(trial);        
        p1 = bpX2pnt(s.window(1) + trialZero, frameRate);
        p2 = bpX2pnt(s.window(2) + trialZero, frameRate);        
        
%         trialZero = zeroTimes2(trial) - pupil.startTime(trial);        
%         p1 = bpX2pnt(s.window(1) + trialZero, frameRate);
%         p2 = bpX2pnt(s.window(2) + trialZero, frameRate);
        
%         p1 = max(1, bpX2pnt(w2(1), frameRate));
%         p2 = min(nFrames, bpX2pnt(w2(2), frameRate));        
        trialData = pupil.(s.pupilField)(trial, p1:p2);
        switch s.method
            case 'mean'
                peak.data(trial) = mean(trialData);
            case 'max'
                peak.data(trial) = max(trialData);
            case 'min'
                peak.data(trial) = min(trialData);                
            otherwise
                
        end
    end