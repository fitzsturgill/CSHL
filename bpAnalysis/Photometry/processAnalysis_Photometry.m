function Photometry = processAnalysis_Photometry(session, varargin)
    % new model of how to pipeline my analysis-  there will only ever be
    % one SessionData per session (it is what Bpod spits out).
    
    % no point of supplying analysisName- from now on just use
    % processSessionsAnalysis_PHotometry whether it is 1 session or
    % multiple
    
%     if isempty(analysisName)
%         analysisName = 'analysis';
%     end
%     if ischar(analysisName)
%         if ~isfield(session, analysisName)
%             session.(analysisName) = struct(); % create a struct with name analysis
%         end
%     end
    SessionData = session.SessionData;
    % this structure will become a field in Session.(analysisName)
    Photometry = struct(...
        'data', [],... % 
        'settings', {varargin}...
        );
    nChannels = 2;  % 2 fluor channels
    blStart = 1; % don't use first X seconds of baseline (to avoid filter onset artifact)
    dFFMode = 'byTrial'; % 'byTrial' or 'bySession', baseline calculated on trial by trial basis or across session
    %default values
    zeroField = '';
    startField = 'PreCsRecording'; % old startField was PreTrialRecording
    trialTypes = unique(SessionData.TrialTypes);
    outcomes = unique(SessionData.TrialOutcome);
    downsample = 305; % factor by which to smooth and decimate
    
    
    % parse input parameter pairs
    counter = 1;
    while counter+1 <= length(varargin) 
        prop = varargin{counter};
        val = varargin{counter+1};
        switch prop
            case 'zeroField'
                zeroField = val;
            case 'startField'
                startField = val;
            case 'trialTypes'
                trialTypes = val;
            case 'outcomes'
                outcomes = val;
            case 'dFFMode'
                dFFMode = val;
            case 'downsample'
                downsample = val;
            otherwise
        end
        counter=counter+2;
    end
    
    if isempty(zeroField)
        disp('*** Error in processAnalysis_Photometry: zeroField must be specified ***');
        return
    end
    
    % make sure photometry data has been demodulated
    if ~isfield(SessionData, 'demod')
        disp('*** First run demodulateSession ***');
        % SessionData = demodulateSession(SessionData);        
    end        

    try
        sampleRate = SessionData.Settings.sample_rate;
    catch
        sampleRate = 6100;
    end    
    sampleRate = sampleRate / downsample;
%     if ~isinteger(sampleRate)
%         disp('error in processAnalysis_Photometry: pick an evenly divisible downsample factor');
%         return
%     end
    trialTypes = unique(SessionData.TrialTypes);
    outcomes = unique(SessionData.TrialOutcome);
    
    % initialize data fields
    data = struct(...
        'data', [],...
        'raw', [],...
        'blF', [],...
        'mod', [],...
        'avg', [],...
        'x', [],...
        'n', [],...
        'std', [],...
        'sem', [],...
        'trials', []...
        );
    Photometry.data = repmat(data, max(trialTypes), max(outcomes));
    
    for i = 1:length(trialTypes)
        type = trialTypes(i);
        for j = 1:length(outcomes)
            outcome = outcomes(j);
            trials = bpFilterTrials(SessionData, type, outcome);
            if isempty(trials)
                continue
            end
            Photometry.data(type, outcome).trials = trials';                    
            for fCh=1:nChannels
        % initialize, assuming that data is of consistent size matching
        % that on the first trial
                originalSamples = max(cellfun(@length, SessionData.demod(:,1))); % just check first channel
                newSamples = ceil(originalSamples/downsample);                
%                 newSamples = ceil(size(SessionData.demod{trials(1), fCh}, 1)/downsample);
                allData = NaN(length(trials), newSamples);
                %modData = raw data that has not been demodulated                
                modData = NaN(length(trials), newSamples); 
                zeroTime = round(SessionData.RawEvents.Trial{trials(1)}.States.(zeroField)(1,1) - SessionData.RawEvents.Trial{trials(1)}.States.(startField)(1,1));
%                 zeroTime = 2;
                avgX = linspace(-zeroTime, ((size(allData, 2) - 1) / sampleRate) - zeroTime, size(allData, 2));
                for counter = 1:length(trials)
                    trial = trials(counter);
                    trialData = SessionData.demod{trial, fCh}'; % convert to row vector
                    %% in case nidaq acquisition ended early for some reason, pad with NaNs
                    if size(trialData, 2) < originalSamples
                        decSamplesShort = newSamples - ceil(size(trialData, 2)/downsample);
                        allData(counter, :) = [decimate(trialData, downsample) NaN(size(trialData, 1), decSamplesShort)]; % pad with NaNs
                        modData(counter, :) = [decimate(SessionData.NidaqData{trial, 1}( :,fCh)', downsample) NaN(size(trialData, 1), decSamplesShort)];                                            
%                         allData(counter, :) = decimate([trialData NaN(size(trialData, 1), samplesShort)], downsample); % pad with NaNs
%                         modData(counter, :) = decimate([SessionData.NidaqData{trial, 1}( :,fCh)' NaN(size(trialData, 1), samplesShort)], downsample);                                                                    
                    else
                        allData(counter, :) = decimate(trialData(1, 1:originalSamples), downsample);
                        modData(counter, :) = decimate(SessionData.NidaqData{trial, 1}( 1:originalSamples,fCh)', downsample);                        
                    end
                end

%             convert to deltaF/F
                blStartP = bpX2pnt(blStart, sampleRate);
                blEndP = bpX2pnt(zeroTime, sampleRate); 
                switch dFFMode
                    case 'byTrial'
                        blF = nanmean(allData(:, blStartP:blEndP), 2); % take mean across time, not trials
                        blF = repmat(blF, 1, size(allData, 2));
                    case 'bySession' % not yet tested
                        blF = nanmean(nanmean(allData(:, blStartP:blEndP), 2)); % scalar mean across time and trials
                        blF = zeros(size(allData)) + blF;                        
                    otherwise
                end
                % rawData not converted to deltaF/F
                rawData = allData;
                allData = (allData - blF) ./ blF;
                avg = nanmean(allData);
                sem = std(allData, 'omitnan') ./ sqrt(sum(~isnan(allData), 1));            
                Photometry.data(type, outcome).data(:,:,fCh) = allData;
                Photometry.data(type, outcome).raw(:,:,fCh) = rawData;       
                Photometry.data(type, outcome).blF(:,fCh) = mean(blF, 2);                    
                Photometry.data(type, outcome).mod(:,:,fCh) = modData;                                
                Photometry.data(type, outcome).avg(:,fCh) = avg';
                Photometry.data(type, outcome).sem(:,fCh) = sem';
                Photometry.data(type, outcome).std(:,fCh) = std(allData, 'omitnan')';
                Photometry.data(type, outcome).n(:,fCh) = sum(~isnan(allData), 1)';
                Photometry.data(type, outcome).x(:,fCh) = avgX';
            end    
        end
    end
    

                
                
                