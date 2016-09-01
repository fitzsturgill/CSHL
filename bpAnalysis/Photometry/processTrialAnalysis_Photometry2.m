function Photometry = processTrialAnalysis_Photometry2(sessions, varargin)
% exemplar for new trial analysis functions
    

    %% optional parameters, first set defaults
    defaults = {...
        'channels', 1;... % 8/28/2016- changed channels default from [] to 1
        'baseline', [1 3];... % 1 - 3 second into recording
        'dFFMode', 'byTrial';...
        'zeroField', 'Us';...
        'startField', 'PreCsRecording';... % TO DO: PROVIDE AUTOMATICALLY BY BPOD NIDAQ CODE 
        'downsample', 305;...
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (store data in cell array)
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    % find total number of trials across selected sessions and size of
    % nidaq data
    scounter = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
    end
    totalTrials = sum(scounter);    
    if s.uniformOutput % not fully implemented
        originalSamples = max(cellfun(@(x) size(x,1), sessions(1).SessionData.NidaqData(:,1)));
        newSamples = ceil(originalSamples/s.downsample);
        try
            sampleRate = sessions(1).SessionData.TrialSettings(1).nidaq.sample_rate;
        catch
            sampleRate = 6100; % very early sessions don't have sample rate in settings
        end
        if rem(sampleRate / s.downsample)
            error('downsample must be a factor of sampleRate');
        end
        sampleRate = sampleRate / s.downsample; % downsample must be a factor of sampleRate      
        zeroTime = sessions(1).SessionData.RawEvents.Trial{1}.States.(s.zeroField)(1) - ...
            sessions(1).SessionData.RawEvents.Trial{1}.States.(s.startField)(1);
        xData = linspace(-zeroTime, ((newSamples - 1) / sampleRate) - zeroTime, newSamples);    
    else
        originalSamples = [];
        newSamples = [];
    end
     

%% Initialize
    Photometry = struct(...
        'data', [],...          % length = number of Channels
        'settings', s,...
        'sampleRate', sampleRate,... % downsampled sample rate
        'startTime', NaN(totalTrials, 1),... % what if a trial didn't have photomtery...., make this NaN initially therefore
        'xData', xData... % you don't have to use xData, you can also use startTime for more flexible alignment to trial events
    );
    if s.uniformOutput
        data = struct(...
            'dFF', NaN(totalTrials, newSamples),... % deltaF/F
            'raw', NaN(totalTrials, newSamples),... %
            'blF', NaN(totalTrials, 1),...
            'ch', []...
            );    
    else
        data = struct(...        
            'dFF', {},... % deltaF/F
            'raw', {},... %
            'blF', NaN(totalTrials, 1),...
            'ch', []...
            );    
    end
    Photometry.data = repmat(data, length(s.channels), 1); % length = # channels                
    h = waitbar(0, 'Processing Photometry');    
    
    tcounter = 1;    
    for si = 1:length(sessions);
        SessionData = sessions(si).SessionData;
        if ~isfield(SessionData, 'demod');
            SessionData = demodulateSession(SessionData); % don't necessarily want to save these back to sessions because that'd eat up memory
        end
        startTimes = cellfun(@(x) x.States.(s.startField)(1), sessions(si).SessionData.RawEvents.Trial); % take the beginning time stamp for the startField-specified Bpod state
        nTrials = SessionData.nTrials;
        allData = NaN(nTrials, newSamples);   
%         modData = NaN(nTrials, newSamples);  % raw data that has not been demodulated   


        for i = 1:length(s.channels)
            fCh = s.channels(i);
            for trial = 1:nTrials
                trialData = SessionData.demod{trial, fCh}'; % convert to row vector
                %% in case nidaq acquisition ended early for some reason, pad with NaNs, this should be fixed as of 8/2016
                if size(trialData, 2) < originalSamples
                    samplesShort = originalSamples - size(trialData, 2);
                    thePad = NaN(size(trialData, 1), samplesShort);
                    allData(trial, :) = decimate([trialData thePad], s.downsample); % pad with NaNs
%                     modData(trial, :) = decimate([SessionData.NidaqData{trial, 1}( :,fCh)' thePad], s.downsample); 
                    disp(['*** short samples on trial ' num2str(trial) ' ***']);
                else
                    allData(trial, :) = decimate(trialData, s.downsample);
%                     modData(trial, :) = decimate(SessionData.NidaqData{trial, 1}(:,fCh)', s.downsample);                        
                end
                                 
            end

            % convert to deltaF/F
            blStartP = bpX2pnt(s.baseline(1), sampleRate);
            blEndP = bpX2pnt(s.baseline(2), sampleRate); 
            switch s.dFFMode
                case 'byTrial'
                    blF = nanmean(allData(:, blStartP:blEndP), 2); % take mean across time, not trials
                    blF = repmat(blF, 1, size(allData, 2));
                case 'bySession' % not yet tested
                    blF = nanmean(nanmean(allData(:, blStartP:blEndP), 2)); % scalar mean across time and trials
                    blF = zeros(size(allData)) + blF;                        
                otherwise
            end      
            Photometry.data(fCh).dFF(tcounter:tcounter+nTrials - 1, :) = (allData - blF) ./ blF;  
            Photometry.data(fCh).raw(tcounter:tcounter+nTrials - 1, :) = allData;                  
            Photometry.data(fCh).blF(tcounter:tcounter+nTrials - 1, 1) = mean(blF, 2);
            Photometry.data(fCh).ch = fCh;
        end
        Photometry.startTime(tcounter:tcounter+nTrials - 1) = startTimes';
        tcounter = tcounter + nTrials;
        waitbar(si/length(sessions));
    end
    close(h);
    
%% Old version      
%     function Photometry = processTrialAnalysis_Photometry(session, varargin)
% 
%     
%     SessionData = session.SessionData;
%     
%     %% optional parameters, first set defaults
%     
%     defaults = {...
%         'nChannels', 2, '';...
%         'blStart', 1, '';...
%         'dFFMode', 'byTrial', '';...
%         'zeroField', 'DeliverStimulus', '';...
%         'startField', 'PreCsRecording','';...
%         'downsample', 100, '';...
%         };
%     settings = bpParseSettings(defaults, varargin); % combine default and passed (via varargin) parameter settings
%     % in future it may be better to use this structure directly.
%     % settings.nChannels.   also the structure you can add to if you like.
%     %e.g. settings.nChannels = SessionData.TrialSettings(1).nidaq.nChannels
%     nChannels = settings.('nChannels');  % 2 fluor channels
%     blStart = settings.('blStart'); % don't use first X seconds of baseline (to avoid filter onset artifact)
%     dFFMode = settings.('dFFMode'); % 'byTrial' or 'bySession', baseline calculated on trial by trial basis or across session
%     %default values
%     zeroField = 'DeliverStimulus';
%     startField = 'PreCsRecording'; % old startField was PreTrialRecording
%     endField = ''; % not implemented
%     downsample = 100; % factor by which to smooth and decimate
%     
%     
%     
%     Photometry = struct(...
%         'data', [],...
%         'settings', settings...
%         );
%     try
%         sampleRate = SessionData.TrialSettings(1).nidaq.sample_rate;
%     catch
%         sampleRate = 6100; % very early sessions don't have sample rate in settings
%     end
%     sampleRate = sampleRate / downsample;
%     
%     data = struct(...
%         'dFF', [],...
%         'raw', [],...
%         'blF', [],...
%         'mod', [],...
%         'x', [],...
%         'n', []...
%         );    
%     Photometry.data = repmat(data, nChannels, 1);
%     for fCh=1:nChannels
%     % initialize, assuming that data is of consistent size matching
%     % that on the first trial
%             originalSamples = size(SessionData.demod{1, fCh}, 1);
%             newSamples = ceil(size(SessionData.demod{1, fCh}, 1)/downsample);
%             allData = NaN(SessionData.nTrials, newSamples);
%             %modData = raw data that has not been demodulated                
%             modData = NaN(SessionData.nTrials, newSamples); 
%             zeroTime = round(SessionData.RawEvents.Trial{1}.States.(zeroField)(1,1) - SessionData.RawEvents.Trial{1}.States.(startField)(1,1));
% %                 zeroTime = 2;
%             x = repmat(linspace(-zeroTime, ((size(allData, 2) - 1) / sampleRate) - zeroTime, size(allData, 2)), SessionData.nTrials, 1);
%         for counter = 1:SessionData.nTrials
%             trial = counter; % copying old code, variable trial is a holdover
%             trialData = SessionData.demod{trial, fCh}'; % convert to row vector
%             %% in case nidaq acquisition ended early for some reason, pad with NaNs
%             if size(trialData, 2) < originalSamples
%                 samplesShort = size(allData, 2) - size(trialData, 2);
%                 allData(counter, :) = decimate([trialData NaN(size(trialData, 1), samplesShort)], downsample); % pad with NaNs
%                 modData(counter, :) = decimate([SessionData.NidaqData{trial, 1}( :,fCh)' NaN(size(trialData, 1), samplesShort)], downsample); 
%                 disp(['shortsamples on trial' num2str(trial)]);
%             else
%                 allData(counter, :) = decimate(trialData(1, 1:originalSamples), downsample);
%                 modData(counter, :) = decimate(SessionData.NidaqData{trial, 1}( 1:originalSamples,fCh)', downsample);                        
%             end
%         end
%         
%         %       convert to deltaF/F
%         blStartP = bpX2pnt(blStart, sampleRate);
%         blEndP = bpX2pnt(zeroTime, sampleRate); 
%         switch dFFMode
%             case 'byTrial'
%                 blF = nanmean(allData(:, blStartP:blEndP), 2); % take mean across time, not trials
%                 blF = repmat(blF, 1, size(allData, 2));
%             case 'bySession' % not yet tested
%                 blF = nanmean(nanmean(allData(:, blStartP:blEndP), 2)); % scalar mean across time and trials
%                 blF = zeros(size(allData)) + blF;                        
%             otherwise
%         end
%         % rawData not converted to deltaF/F
%         rawData = allData;
%         allData = (allData - blF) ./ blF;          
%         
%         Photometry.data(fCh).dFF = (allData - blF) ./ blF;  
%         Photometry.data(fCh).raw = allData;                  
%         Photometry.data(fCh).blF = mean(blF, 2);
%         Photometry.data(fCh).mod = modData;
%         Photometry.data(fCh).x = x;
%     end
%     
%             
%     
%     