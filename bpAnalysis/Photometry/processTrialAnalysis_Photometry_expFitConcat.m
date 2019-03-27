function Photometry = processTrialAnalysis_Photometry_expFitConcat(sessions, varargin)
% Fitz Sturgill 2018
% This function subtracts away an exponential bleaching trend by
% concatenating all trials together.  Assumes that bleaching occurs approximately without
% replenishment of unbleached fluorophore from a pool outside the
% photometry "field of view".  This assumtion is supported by the
% observation that fluorecense tends not to recover across successive days
% of photometry sessions.

% unlike processTrialAnalysis_Photometry2, this function cuts away
% the options for other types of bleaching correction and baselining for simplicity    

    %% optional parameters, first set defaults
    defaults = {...
        'channels', 1;... % 8/28/2016- changed channels default from [] to 1
        'refChannels', [];...
        'baseline', [1 4];... % 1 - 3 second into recordinge
        'zeroField', 'Us';...
        'startField', 'PreCsRecording';... % TO DO: PROVIDE AUTOMATICALLY BY BPOD NIDAQ CODE 
        'downsample', 305;...
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (store data in cell array)
        'forceAmp', 0;... % % force demodulation even if the refChannel LED is off (i.e. it's amplitude = 0)
        'ACfilter', 0;...
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    %% if not specified per channel (as a cell array), use identical values across channels

    if ~iscell(s.baseline)
        s.baseline = repmat({s.baseline}, 1, max(s.channels));
    end
    
    if isempty(s.refChannels)
        s.refChannels = s.channels; % use same reference by default
    end   
    
    % find total number of trials across selected sessions and size of
    % nidaq data
    scounter = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
    end
    totalTrials = sum(scounter);    
    
    if s.uniformOutput % different sized trials are padded with NaNs, aligned to photometry start
        % determine maximum number of samples
        maxDuration = 0;
        for scounter = 1:length(sessions)
            for counter = 1:length(sessions(scounter).SessionData.TrialSettings)
                maxDuration = max(maxDuration, sessions(scounter).SessionData.TrialSettings(counter).nidaq.duration);
            end
        end
        try
            sampleRate = sessions(1).SessionData.TrialSettings(1).nidaq.sample_rate;
        catch
            sampleRate = 6100; % very early sessions don't have sample rate in settings
        end
        originalSamples = maxDuration * sessions(1).SessionData.TrialSettings(1).nidaq.sample_rate;                            
        newSamples = ceil(originalSamples/s.downsample);
        if rem(sampleRate, s.downsample)
            error('downsample must be a factor of sampleRate');
        end
        % update sampleRate because you no longer need non-decimated sample rate
        sampleRate = sampleRate / s.downsample; % downsample must be a factor of sampleRate      
        zeroTime = sessions(1).SessionData.RawEvents.Trial{1}.States.(s.zeroField)(1) - ...
            sessions(1).SessionData.RawEvents.Trial{1}.States.(s.startField)(1);
        xData = linspace(-zeroTime, ((newSamples - 1) / sampleRate) - zeroTime, newSamples);    
    else % unfinished, different sized trials are stored in cell arrays
        originalSamples = [];
        newSamples = [];
    end
     

%% Initialize
    Photometry = struct(...
        'data', [],...          % length = number of Channels
        'settings', s,...
        'bleachFit', [],... % length
        'sampleRate', sampleRate,... % downsampled sample rate
        'startTime', NaN(totalTrials, 1),... % what if a trial didn't have photometry...., make this NaN initially therefore
        'xData', xData... % you don't have to use xData, you can also use startTime for more flexible alignment to trial events
    );
    if s.uniformOutput % different sized trials are padded with NaNs, aligned to photometry start
        data = struct(...
            'dFF', NaN(totalTrials, newSamples),... % deltaF/F
            'dF', NaN(totalTrials, newSamples),... % deltaF            
            'ZS', NaN(totalTrials, newSamples),... % deltaF/F            
            'raw', NaN(totalTrials, newSamples),... %
            'blF_fit', NaN(totalTrials, newSamples),...
            'ch', []...
            );    
    else % unfinished, different sized trials are stored in cell arrays
        data = struct(...        
            'dFF', {},... % deltaF/F
            'dF', {},...
            'ZS', {},... % deltaF/F            
            'raw', {},... %
            'blF_fit', {},...
            'ch', []...
            );    
    end
    
    bleachFit = struct(... % outputs from fit function in curve fitting toolbox
        'fitobject_session', [],...
        'gof_session', [],...
        'output_session', []...
        );

    Photometry.data = repmat(data, length(s.channels), 1); % length = # channels  
    Photometry.bleachFit = repmat(bleachFit, length(sessions), max(s.channels));
    h = waitbar(0, 'Processing Photometry');    

    
    tcounter = 1;    
    for si = 1:length(sessions)
        SessionData = sessions(si).SessionData;
        if ~isfield(SessionData, 'demod')
            SessionData = demodulateSession(SessionData, 'channels', s.channels, 'refChannels', s.refChannels, 'forceAmp', s.forceAmp, 'ACfilter', s.ACfilter); % don't necessarily want to save these back to sessions because that'd eat up memory
        end
        startTimes = cellfun(@(x) x.States.(s.startField)(1), sessions(si).SessionData.RawEvents.Trial); % take the beginning time stamp for the startField-specified Bpod state
        nTrials = SessionData.nTrials;
        allData = NaN(nTrials, newSamples);    

        for chCounter = 1:length(s.channels)
            fCh = s.channels(chCounter);
            for trial = 1:nTrials
                trialData = SessionData.demod{trial, fCh}'; % convert to row vector
                %% in case nidaq acquisition ended early for some reason, pad with NaNs, this should be fixed as of 8/2016
                downData = decimate(trialData, s.downsample);
                if length(downData) < newSamples
                    downData = [downData NaN(1, newSamples - length(downData))];
                elseif length(downData) > newSamples
                    downData = downData(1:newSamples);
                end
                allData(trial, :) = downData;
            end

            % baseline fluor.
            blStartP = bpX2pnt(s.baseline{chCounter}(1), sampleRate);
            blEndP = bpX2pnt(s.baseline{chCounter}(2), sampleRate); 
            
            
            % concatenate all data together
            blF_fit = allData';
            blF_fit = blF_fit(:);
%             fo = fitoptions('Method', 'NonlinearLeastSquares',...
%                 'Upper', [Inf range(blF_fit) 0 range(blF_fit) 0],...
%                 'Lower', [0 0 -1 0 -1],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
%                 'StartPoint', [min(blF_fit) range(blF_fit)/2 -5 range(blF_fit)/2 -100]...
%                 );
%             model = 'a + b*exp(c*x) + d*exp(e*x)';
            fo = fitoptions('Method', 'NonlinearLeastSquares',...
                'Upper', [Inf range(blF_fit) 0],...
                'Lower', [0 0 -1],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
                'StartPoint', [min(blF_fit) range(blF_fit)/2 -5]...
                );
            model = 'a + b*exp(c*x)';
            ft = fittype(model, 'options', fo);
            if any(isnan(blF_fit))
                blF_fit = inpaint_nans(blF_fit);
            end
            x = (0:length(blF_fit)-1)';
            [fitobject, gof, output] = ... % new
                fit(x, blF_fit, ft, fo);
            Photometry.bleachFit(si, fCh).fitobject_session = fitobject;
            Photometry.bleachFit(si, fCh).gof_session = gof;
            Photometry.bleachFit(si, fCh).output_session = output;

%             blF_fit = fitobject.a + fitobject.b * exp(fitobject.c * x) + fitobject.d * exp(fitobject.e * x);   
            blF_fit = fitobject.a + fitobject.b * exp(fitobject.c * x);
            
            % subtract away the bleaching trend
            blF_fit = reshape(blF_fit, size(allData, 2), size(allData, 1))';
            dF = allData - blF_fit;
 
            Photometry.data(fCh).dFF(tcounter:tcounter+nTrials - 1, :) = dF ./ blF_fit; 
            Photometry.data(fCh).dF(tcounter:tcounter+nTrials - 1, :) = dF;
            sd = nanmean(nanstd(dF(:,blStartP:blEndP), 0, 2));
            Photometry.data(fCh).ZS(tcounter:tcounter+nTrials - 1, :) = dF / sd; % z-scored
            Photometry.data(fCh).raw(tcounter:tcounter+nTrials - 1, :) = allData;                   
            Photometry.data(fCh).blF_fit(tcounter:tcounter+nTrials - 1, :) = blF_fit;
            Photometry.data(fCh).ch = fCh;
        end
        Photometry.startTime(tcounter:tcounter+nTrials - 1) = startTimes';
        tcounter = tcounter + nTrials;
        waitbar(si/length(sessions));
    end
    close(h);
    

