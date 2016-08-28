function Photometry = processTrialAnalysis_Photometry2(sessions, Photometry, varargin)
% exemplar for new trial analysis functions
    

    %% optional parameters, first set defaults
    defaults = {...
        'channels', [];...
        'blStart', 1;...
        'dFFMode', 'byTrial';...
        'zeroField', 'DeliverStimulus';...
        'startField', 'PreCsRecording';... % TO DO: PROVIDE AUTOMATICALLY BY BPOD NIDAQ CODE 
        'downsample', 500;...
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (store data in cell array)
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings

%% Initialize
    for si = 1:length(sessions);
        SessionData = session(si).SessionData;
        try
            sampleRate = SessionData.TrialSettings(1).nidaq.sample_rate;
        catch
            sampleRate = 6100; % very early sessions don't have sample rate in settings
        end
        sampleRate = sampleRate / s.downsample; % I'm assuming downsample is a factor of sampleRate
        startTimes = cellfun(@(x) x.States.(s.startField)(1), session.SessionData.RawEvents.Trial); % take the beginning time stamp for the startField-specified Bpod state

        Photometry = struct(...
            'data', [],...          % length = number of Channels
            'settings', s,...
            'sampleRate', sampleRate,... % downsampled sample rate
            'startTime',startTimes...
        );

        data = struct(...
            'dFF', [],... %deltaF/F
            'raw', [],... %
            'blF', [],...
            'mod', []...
            );    
        Photometry.data = repmat(data, s.nChannels, 1); % length = # channels

        %% initialize
        originalSamples = max(cellfun(@length, SessionData.demod(:,1))); % just check first channel
        newSamples = ceil(originalSamples/s.downsample);
        allData = NaN(SessionData.nTrials, newSamples);   
        modData = NaN(SessionData.nTrials, newSamples);  % raw data that has not been demodulated         

        for fCh=1:s.nChannels
            if s.constantLength
                % baselineEndTime is relative to start of photometry recording
                baselineEndTime = round(SessionData.RawEvents.Trial{1}.States.(s.zeroField)(1,1) - SessionData.RawEvents.Trial{1}.States.(s.startField)(1,1));
            else
                disp('*** variable trial length not yet implemented ***');
            end
            for trial = 1:SessionData.nTrials
                trialData = SessionData.demod{trial, fCh}'; % convert to row vector
                %% in case nidaq acquisition ended early for some reason, pad with NaNs
                %   WHY DOES THIS STILL HAPPEN????
                if size(trialData, 2) < originalSamples
                    samplesShort = originalSamples - size(trialData, 2);
                    thePad = NaN(size(trialData, 1), samplesShort);
                    allData(trial, :) = decimate([trialData thePad], s.downsample); % pad with NaNs
                    modData(trial, :) = decimate([SessionData.NidaqData{trial, 1}( :,fCh)' thePad], s.downsample); 
                    disp(['*** short samples on trial ' num2str(trial) ' ***']);
                else
                    allData(trial, :) = decimate(trialData, s.downsample);
                    modData(trial, :) = decimate(SessionData.NidaqData{trial, 1}(:,fCh)', s.downsample);                        
                end
            end

            %       convert to deltaF/F
            blStartP = bpX2pnt(s.blStart, sampleRate);
            blEndP = bpX2pnt(baselineEndTime, sampleRate); 
            switch s.dFFMode
                case 'byTrial'
                    blF = nanmean(allData(:, blStartP:blEndP), 2); % take mean across time, not trials
                    blF = repmat(blF, 1, size(allData, 2));
                case 'bySession' % not yet tested
                    blF = nanmean(nanmean(allData(:, blStartP:blEndP), 2)); % scalar mean across time and trials
                    blF = zeros(size(allData)) + blF;                        
                otherwise
            end      
            Photometry.data(fCh).dFF = (allData - blF) ./ blF;  
            Photometry.data(fCh).raw = allData;                  
            Photometry.data(fCh).blF = mean(blF, 2);
            Photometry.data(fCh).mod = modData;
        end
    end
    
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