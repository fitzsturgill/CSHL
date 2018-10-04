function peak = bpCalcPeak_Pupil(pupil, varargin)
% pupil- output of addPupilometryToTE
% Can handle differently sized whisk data per trial (nonUniformOutput mode)

% see also bpCalcPeak_dFF

    defaults = {...
        'method', 'mean';... % mean | min | max | percentile (in which case percentile must be specified)
        'pupilField', 'pupDiameterNorm';...
        'zeroTimes', [];... % 1) empty- inferred from whisk xData 2) cell array or 3) vector containing zero times relative to Bpod trial start
        'window', [];... % REQUIRED, averaging window around alignment point
        'referenceFromEnd', 0;... % relevent when zeroTimes are supplied as a cell array (e.g. if you want to align to the beginning or end of a bpod state)        
        'percentile', [];...
        };
    
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    assert(~isempty(s.window), 'Must provide window for whisk peak calculation');

    if strcmp(s.method, 'percentile')
        assert(~isempty(s.percentile), 'percentile must be specified');
    end
    nTrials = size(pupil.startTime, 1);

    if isnumeric(pupil.(s.pupilField))
        uniformOutput = 1;
    elseif iscell(pupil.(s.pupilField))
        uniformOutput = 0;
    else
        error('pupil data of improper type');
    end
    peak = struct(...
        'data', zeros(nTrials, 1),...
        'settings', s...
    );

    if isempty(s.zeroTimes)
        s.zeroTimes = pupil.startTime - pupil.xData(1);
    elseif iscell(s.zeroTimes)
        if ~s.referenceFromEnd % first point in each cell array element
            zeroTimes = cellfun(@(x) x(1), s.zeroTimes); % returns vector
        else % last point in each cell array element
            zeroTimes = cellfun(@(x) x(end), zeroTimes); % returns vector
        end
    else
        zeroTimes = s.zeroTimes; 
    end
    
    zeroTimes = zeroTimes(:); 
    
    if size(s.window, 1) == 1
        s.window = repmat(s.window, nTrials, 1);
    end
    
    zeroTimes = zeroTimes - pupil.startTime; % redefine zero times relative to pupil start    
    
    for trial = 1:nTrials
        sampleRate = pupil.frameRate(trial);    
        if uniformOutput
            trialData = pupil.(s.pupilField)(trial, :);
        else
            trialData = pupil.(s.pupilField){trial};
        end
        nSamples = length(trialData);
        acqDuration = nSamples / sampleRate;
        zeroTime = zeroTimes(trial);
        if isnan(zeroTime) || (zeroTime < 0) || (zeroTime > acqDuration)|| any(isnan(s.window(trial, :)))
            continue % continue if zero time occurs outside photometry acquisition or if windows or zeros contain NaNs
        end        
        p1 = bpX2pnt(zeroTime + s.window(trial, 1), sampleRate);
        p2 = min(bpX2pnt(zeroTime + s.window(trial, 2), sampleRate), nSamples);        

        switch s.method
            case 'mean'
                peak.data(trial) = mean(trialData(p1:p2));
            case 'max'
                peak.data(trial) = max(trialData(p1:p2));
            case 'min'
                peak.data(trial) = min(trialData(p1:p2));                
            case 'percentile'
                peak.data(trial) = percentile(trialData(p1:p2), s.percentile);
            otherwise
        end
    end