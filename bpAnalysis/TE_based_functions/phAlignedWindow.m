function [data, xData] = phAlignedWindow(TE, trials, ch, varargin)

    defaults = {...
        'PhotometryField', 'Photometry';...
        'FluorDataField', 'dFF';...
        'window', [];... % averaging window around alignment point
        'zeroTimes', [];... % empty- inferred from Photometry xData OR length TE-nTrials cell array or vector containing zero times relative to Bpod trial start
        'referenceFromEnd', 0;... % relevent when zeroTimes are supplied as a cell array (e.g. if you want to align to the beginning or end of a bpod state)        
        };    
    [s, ~] = parse_args(defaults, varargin{:});


    if ~isfield(TE, s.PhotometryField)
        error([Photometry ' field does not exist']);
    else
        Photometry = TE.(s.PhotometryField);
    end
    
    Fs = 1/Photometry.sampleRate;

    if islogical(trials)
        trials = find(trials);
    end
    nTrials = length(trials);

    
%%
    if isempty(s.zeroTimes)
        zeroTimes = valuecrossing(1:length(TE.Photometry.xData), TE.Photometry.xData, 0); % inferred from xData, historical usage...
        if isempty(s.window)
            s.window = TE.Photometry.xData([1 end]);
        end
    elseif iscell(s.zeroTimes)
        if ~s.referenceFromEnd
            zeroTimes = cellfun(@(x) x(1), s.zeroTimes); % returns vector
        else
            zeroTimes = cellfun(@(x) x(end), zeroTimes); % returns vector
        end
        assert(~isempty(s.window), 'user must supply window if zero time(s) are supplied as well');
    else
        zeroTimes = s.zeroTimes; 
        assert(~isempty(s.window), 'user must supply window if zero time(s) are supplied as well');
    end
    
    zeroTimes = zeroTimes(:); 
    if isscalar(zeroTimes)
        zeroTimes = repmat(zeroTimes, nTrials, 1);
    end
        
    if size(s.window, 1) == 1
        s.window = repmat(s.window, nTrials, 1);
    end
    
    zeroTimes = zeroTimes - TE.Photometry.startTime; % redefine zero times relative to photometry start
    

%% determine maximum number of data points preceding and following zero and initalize data array
    if Photometry.settings.uniformOutput
        samplesPerTrial = size(TE.Photometry.data(ch).raw); % how many samples per trial
        paddedSamples = bpX2pnt(max(zeroTimes)) + samplesPerTrial - bpX2pnt(min(zeroTimes)); % for padding array with maximum number of points before and after a trial zero
        paddedZeroPoint = bpX2pnt(max(zeroTimes));
        xData = (0:(samplesPerTrial - 1))/Fs - paddedZeroPoint;
    else
        % FINISH CODING!!!!
        error('finish coding to allow nonUniformOutput');
    end
    
    data = NaN(nTrials, paddedSamples);
    
%% loop over trials to fill in data array
    for counter = 1:nTrials
        trial = trials(counter);
        zeroTime = zeroTimes(trial);
        if Photometry.settings.uniformOutput
            trialData = Photometry.data(ch).(s.FluorDataField)(trial, :);
            nSamples = length(trialData);
            acqDuration = nSamples / Fs;
        else
            % FINISH CODING!!!!
        end
        if (zeroTime < 0) || (zeroTime > acqDuration)
            continue % continue if zero time occurs outside photometry acquisition
        end
        % source data start, end, and zero points
        szero = bpX2pnt(zeroTime, Fs);
        s1 = bpX2pnt(zeroTime + s.window(1), Fs);
        s2 = min(bpX2pnt(zeroTime + s.window(2), Fs), nSamples);
        % destination data start and end points
        d1 = paddedZeroPoint - (szero - s1);
        d2 = paddedZeroPoint + (s2 - szero);
        % add trial data to padded array
        data(counter, d1:d2) = trialData(s1:s2);
    end
        