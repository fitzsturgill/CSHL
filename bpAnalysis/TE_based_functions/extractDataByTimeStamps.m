function varargout = extractDataByTimeStamps(data, startTimes, Fs, TS, window, trials)
% outputs [eventData, timeStamps, trialNumbers];
% extracts time stamp triggered windows of continuous data derived from a
% TE structure containing data from dataField contains following field
% data- continuous data, nTrials x nSamples or nTrials x nSamplesY x
% nSamplesZ (for 3 dimensional data, 2nd dimension size is == to number of
% time points
% startTimes- startTimes of continuous data relative to trial start (scalar or nTrials x 1)
% Fs- sample rate of continuous data
% TS, time stamps, cell array nTrials x 1, assumes that time stamps are in
% the first column of cell array contents for a given trial (so that you
% can use output of bpAddStateAsTrialEvent directly)
% window, time window around time stamp to extract
    if nargin < 6
        trials = 1:length(TS);
    end
    



% TE.Photometry.data(1).ZS
% TE.Wheel.data.V
% fitz = cellfun(@(x) x(:,1), TE.Reward, 'UniformOutput', 0);
    validTrials = cellfun(@(x) ~isnan(x(:,1)), TS, 'UniformOutput', 0); % use time stamps in first column
    
    validTrials = cellfun(@(x) sum(x), validTrials);    
    if islogical(trials)
        trials = find(trials);
    end
    nWindows = sum(validTrials(trials));
    samplesPerWindow = floor((window(2) - window(1)) * Fs);
    eventData = NaN(nWindows, samplesPerWindow, size(data, 3)); % singleton 3rd dimension in case of 2d data
    
    dT = 1/Fs;    
    eventCounter = 1;
    validTrials = intersect(find(validTrials), trials);
    validTrials = validTrials(:)';
    timeStamps = zeros(nWindows, 1);
    trialNumbers = zeros(nWindows, 1);
    for trial = validTrials
        startTime = startTimes(trial);
        if iscell(data)
            warning('Take a look at code below, do I need to use shiftdim, why do I use a colon in the 3rd dimension?');
            trialData = data{trial}(:,:,:);
        else
            trialData = data(trial, :, :);
            if size(trialData, 1) > 1
                trialData = shiftdim(trialData, -1); % add singleton leading dimension in case you are dealing with 3d data, i.e. 2d data for a given trial
            end
        end
        nPoints = size(trialData, 2);
        stamps = (TS{trial}(:,1)); % use first column of time stamps (in case you are using output of bpAddStateAsTrialEvent);
        stamps = stamps(:)'; % ensure row
        for stamp = stamps
            rs = stamp - startTime; % rs, Relative time Stamp
            trialWindow = window + rs;
            startP = round(1 + trialWindow(1)/dT); % see bpX2Pnt
            sourcePoints = startP:startP + (samplesPerWindow - 1);
            destPoints = find(sourcePoints >= 1 & sourcePoints <= nPoints);
            eventData(eventCounter, destPoints, :) = trialData(1, sourcePoints(destPoints), :);
            timeStamps(eventCounter) = stamp;
            trialNumbers(eventCounter) = trial;
            eventCounter = eventCounter + 1;
        end
    end
            
    if nargout > 0
        varargout{1} = eventData;
    end
    
    if nargout > 1
        varargout{2} = timeStamps;
    end
    
    if nargout > 2
        varargout{3} = trialNumbers;
    end
    
    if nargout > 3
        timeVector = 0:dT:(samplesPerWindow - 1) * dT;
        varargout{4} = timeVector + window(1);
    end
    
    
    