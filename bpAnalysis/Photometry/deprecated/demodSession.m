function [avg, avgX, avgSEM] = demodSession(SessionData, type, outcome, fCh, zeroField, ax, linespec, sampleRate)
    % FS NOTE-  1/13/16,   this is a rough draft, currenty assumes that all
    % fluo data is same length and alignment and doesn't establish zero
    % time


    
% 
    if nargin < 8
        try
            sampleRate = SessionData.Settings.sample_rate;
        catch
            sampleRate = 6100;
        end
    end
    
    if nargin < 5
        zeroField = '';
    end
    
    if nargin < 6
        ax = []; % 1/13/16- haven't completed/decided optional parameters
    end
    
    modF = SessionData.TrialSettings(1,1).(['LED' num2str(fCh) '_f']); % this could change (I should store modF in settings not trialsettings)
    tBaseline = SessionData.TrialSettings(1,1).PreTrialRecording; % use pretrial recording period for normalization/zscoring
    
    
    
    trials = bpFilterTrials(SessionData, type, outcome);
%     zeroField = 'DeliverStimulus'; % this could change
%     startField = 'ITI';
%     endField = 'PostTrialRecording';

%     dataLim = [-2, 4];
    % initialize NaN matrix and xdata

%     avgX = linspace(dataLim(1), dataLim(2), size(SessionData, 2));
    
    %gather trials and align to zero
    for counter = 1:length(trials)
        trial = trials(counter);

        rawData = SessionData.NidaqData{trial,1}(:,fCh);
        nSamples = size(rawData, 1);
        refData = SessionData.NidaqData{trial,2}(1:nSamples,fCh);
%         [finalData, dmf, dmp] = decode_lockin_fn(rawData, refData, [], modF, sampleRate, 0);
        finalData = phDemod(rawData, refData, sampleRate, modF, 30); % 30Hz lowpass corner freq
        
%         trialData = SessionData.NidaqData{trial}(:, fCh)';
%         trialData = decimate(trialData, dmf); % decimate
%         trialData = smooth(trialData, 50, 'lowess')'; % smooth
%         i1 = round((zeroTime + dataLim(1)) * sample_rate);
        % I'm only dealing with trials that don't have enough data after
        % zerotime, not ones that begin too close to zero time (for now):
%         i2 = min(i1 + size(SessionData, 2) - 1, size(trialData, 2));   
        if counter == 1
            allData = NaN(length(trials), size(finalData, 1));
            avgX = linspace(0, nSamples/sampleRate, nSamples);
            if ~isempty(zeroField) 
                zeroTime = SessionData.RawEvents.Trial{trial}.States.(zeroField)(1,1); % always use first instance and start of zeroField for zeroing (each state has a start and end timestamp
                avgX = avgX - zeroTime;
            end
        end
        finalData = finalData'; % convert to row vector
        allData(counter, :) = finalData;
    end
    
    % convert to deltaF/F
%     baselines = nanmean(SessionData(:, 1: abs(dataLim(1)) * sample_rate), 2);
%     baselines = repmat(baselines, 1, size(SessionData, 2));
%     SessionData = SessionData - baselines;
%     SessionData = SessionData ./ baselines;
%     
    
    % avg
    avg = nanmean(allData);
    avgSEM = std(allData, 'omitnan') ./ sqrt(sum(~isnan(allData), 1));

    % only plot data if given axis
    if ~isempty(ax) 
        boundedline(avgX, avg, avgSEM, linespec, ax, 'alpha');
    end
        
        
        
        
        