function Photometry = processTrialAnalysis_Photometry2(sessions, varargin)
% exemplar for new trial analysis functions
    

    %% optional parameters, first set defaults
    defaults = {...
        'channels', 1;... % 8/28/2016- changed channels default from [] to 1
        'refChannels', [];...
        'baseline', [1 4];... % 1 - 3 second into recording
        'blMode', 'byTrial';... % 'byTrial', 'bySession', 'expFit'  expFit- interpolates baseline from biexponential fit to raw fl baselines across trials
        'dFFMode', 'simple';... % 'simple', 'expFit'   !(now with hard-coded time constant for exponential) expFit- subtracts within-trial exponential bleaching trend using an exponential fit to the trial average baseline period
        'expFitBegin', 0.1;...
        'zeroField', 'Us';...
        'startField', 'PreCsRecording';... % TO DO: PROVIDE AUTOMATICALLY BY BPOD NIDAQ CODE 
        'downsample', 305;...
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (store data in cell array)
        'tau', 3;...
        'forceAmp', 0;... % % force demodulation even if the refChannel LED is off (i.e. it's amplitude = 0)
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    %% if not specified per channel (as a cell array), use identical valies across channels
    if ~iscell(s.blMode)
        s.blMode = repmat({s.blMode}, 1, max(s.channels));
    end
    if ~iscell(s.dFFMode)
        s.dFFMode = repmat({s.dFFMode}, 1, max(s.channels));
    end
    if ~iscell(s.baseline)
        s.baseline = repmat({s.baseline}, 1, max(s.channels));
    end
    
    if ~iscell(s.expFitBegin)
        s.expFitBegin = repmat({s.expFitBegin}, 1, max(s.channels));
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
    if s.uniformOutput % not fully implemented
%         originalSamples = zeros(1,length(sessions));
%         for si = 1:length(sessions);
%             originalSamples(si) = max(cellfun(@(x) size(x,1), sessions(si).SessionData.NidaqData(:,1)));
%         end
%         originalSamples = max(originalSamples);

        try
            sampleRate = sessions(1).SessionData.TrialSettings(1).nidaq.sample_rate;
            originalSamples = sessions(1).SessionData.TrialSettings(1).nidaq.duration * sessions(1).SessionData.TrialSettings(1).nidaq.sample_rate;            
        catch
            sampleRate = 6100; % very early sessions don't have sample rate in settings
            originalSamples = sessions(1).SessionData.TrialSettings(1).nidaq.duration * sampleRate;            
        end
        
        newSamples = ceil(originalSamples/s.downsample);
        if rem(sampleRate, s.downsample)
            error('downsample must be a factor of sampleRate');
        end
        % update sampleRate because you no longer need non-decimated sample
        % rate
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
        'bleachFit', [],... % length
        'sampleRate', sampleRate,... % downsampled sample rate
        'startTime', NaN(totalTrials, 1),... % what if a trial didn't have photometry...., make this NaN initially therefore
        'xData', xData... % you don't have to use xData, you can also use startTime for more flexible alignment to trial events
    );
    if s.uniformOutput
        data = struct(...
            'dFF', NaN(totalTrials, newSamples),... % deltaF/F
            'dF', NaN(totalTrials, newSamples),... % deltaF            
            'ZS', NaN(totalTrials, newSamples),... % deltaF/F            
            'raw', NaN(totalTrials, newSamples),... %
            'blF', NaN(totalTrials, 1),...
            'blF_raw', NaN(totalTrials, 1),...
            'ch', []...
            );    
    else
        data = struct(...        
            'dFF', {},... % deltaF/F
            'dF', {},...
            'ZS', {},... % deltaF/F            
            'raw', {},... %
            'blF', NaN(totalTrials, 1),...
            'blF_raw', NaN(totalTrials, 1),...
            'ch', []...
            );    
    end
    
    bleachFit = struct(... % outputs from fit function in curve fitting toolbox
        'fitobject_session', [],...
        'gof_session', [],...
        'output_session', [],...
        'fitobject_trial', [],...
        'gof_trial', [],...
        'output_trial', [],...
        'trialFit', [],... % doesn't exclude expFitBegin points
        'trialTemplate', [],... 
        'trialTemplateFull', [],... % doesn't exclude expFitBegin points 
        'trialTemplateX', [],... 
        'trialTemplateFullX', [],... % doesn't exclude expFitBegin points        
        'fitX', []...        % doesn't exclude expFitBegin points
        );

    Photometry.data = repmat(data, length(s.channels), 1); % length = # channels  
    Photometry.bleachFit = repmat(bleachFit, length(sessions), max(s.channels));
    h = waitbar(0, 'Processing Photometry');    
    
    % KLUDGE 
    warning('kludgy implementation of tau per channel, will break if you use only channel 2');
    if length(s.channels) > 1 && length(s.tau) == 1
        s.tau = repmat(s.tau, 1, length(s.channels));
    end
    tcounter = 1;    
    for si = 1:length(sessions)
        SessionData = sessions(si).SessionData;
        if ~isfield(SessionData, 'demod')
            SessionData = demodulateSession(SessionData, 'channels', s.channels, 'refChannels', s.refChannels, 'forceAmp', s.forceAmp); % don't necessarily want to save these back to sessions because that'd eat up memory
        end
        startTimes = cellfun(@(x) x.States.(s.startField)(1), sessions(si).SessionData.RawEvents.Trial); % take the beginning time stamp for the startField-specified Bpod state
        nTrials = SessionData.nTrials;
        allData = NaN(nTrials, newSamples);   
%         modData = NaN(nTrials, newSamples);  % raw data that has not been demodulated   

        for fCh = s.channels
            blMode = s.blMode{fCh};
            dFFMode = s.dFFMode{fCh};
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
            blStartP = bpX2pnt(s.baseline{fCh}(1), sampleRate);
            blEndP = bpX2pnt(s.baseline{fCh}(2), sampleRate); 
            expFitStartP = bpX2pnt(s.expFitBegin{fCh}, sampleRate);
            switch blMode
                case 'byTrial'
                    blF = nanmean(allData(:, blStartP:blEndP), 2); % take mean across time, not trials
                    blF_raw = blF;                    
                    blF = repmat(blF, 1, size(allData, 2));
                case 'bySession' % not yet tested
                    blF = nanmean(nanmean(allData(:, blStartP:blEndP), 2)); % SCALAR mean across time and trials
                    blF = zeros(size(allData)) + blF;
                    blF_raw = mean(blF, 2);
                case 'expFit'
                    blF_raw = nanmean(allData(:, expFitStartP:blEndP), 2); % take mean across time, not trials
%                     blF_fit = medfilt1(blF_raw, 3, 'truncate'); % median filter baseline fluorescence
                    blF_fit = blF_raw;
                    fo = fitoptions('Method', 'NonlinearLeastSquares',...
                        'Upper', [Inf range(blF_raw) 0 range(blF_raw) 0],...
                        'Lower', [0 0 -1 0 -1],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
                        'StartPoint', [min(blF_raw) range(blF_raw)/2 -5 range(blF_raw)/2 -100]...
                        );
                    model = 'a + b*exp(c*x) + d*exp(e*x)';
                    ft = fittype(model, 'options', fo);
                    if any(isnan(blF_fit))
                        blF_fit = inpaint_nans(blF_fit);
                    end
                    [fitobject, gof, output] = ...
                        fit((1:size(allData, 1))', blF_fit, ft, fo);
                    Photometry.bleachFit(si, fCh).fitobject_session = fitobject;
                    Photometry.bleachFit(si, fCh).gof_session = gof;
                    Photometry.bleachFit(si, fCh).output_session = output;
                    x = (0:length(blF_fit)-1)';
                    blF = fitobject.a + fitobject.b * exp(fitobject.c * x) + fitobject.d * exp(fitobject.e * x);
                    blF = repmat(blF, 1, size(allData, 2));
                otherwise
            end
%%        commented code below is snippet to show what different coefficients do to an exponential
% note! c in the example below is the inverse of the time constant, a is
% the asymptote at x = Inf, a + b is the y intercept
% x = 1:1000;
% 
% a = 100;
% b = 1;
% c = -0.01;
% 
% ensureFigure('test', 1);
% plot(x, a + b*exp(c * x), 'k'); hold on;
% plot(x, a + b * 2 *exp(c * x), 'b');
% plot(x, a + b * exp(c * 2 * x), 'r');
% plot(x, a * 1.01 + b * exp(c * x), 'y');
%%
            switch dFFMode
                case 'simple'
                    dF = allData - blF;
                case 'expFit'
                    trialMeanY = (nanmean(allData(:, expFitStartP:blEndP), 1))';
                    trialMeanX = ((0:length(trialMeanY) - 1)/sampleRate + expFitStartP/sampleRate)';   
                    trialMeanYFull = (nanmean(allData(:, 1:blEndP), 1))';        
                    x = (0:size(allData, 2) - 1)/sampleRate;
%% biexponential
%                     fo = fitoptions('Method', 'NonlinearLeastSquares',...
%                         'Lower', [0, 0, -Inf, 0 -Inf],...
%                         'Upper', [min(trialMeanY) Inf 0 Inf 0],...
%                         'StartPoint', [min(trialMeanY) * 0.99, 0.01, -0.1, 0.01, -0.1]);                    
%                     ft = fittype('a + b*exp(c*x) + d*exp(e*x)', 'options', fo);
%                     [fitobject, gof, output] = ...
%                         fit(trialMeanX, trialMeanY, ft, fo);
%                     trialFit = fitobject.a + fitobject.b * exp(fitobject.c * x) + fitobject.d * exp(fitobject.e * x);
%% single exponential with fixed time constant
                    tau = s.tau(fCh); % make this an option later, tau = time constant
                    c = 1/tau;
                    fo = fitoptions('Method', 'NonlinearLeastSquares',...
                        'Lower', [0, range(trialMeanY) * 0.25],...%, -Inf],...
                        'Upper', [min(trialMeanY) Inf],... 0],...
                        'StartPoint', [min(trialMeanY) * 0.99, 0.01]);%, -0.3,]);
%                     fo = fitoptions('Method', 'NonlinearLeastSquares',...
%                         'Lower', [0, 0, -1],...
%                         'Upper', [min(trialMeanY) Inf 0],...
%                         'StartPoint', [min(trialMeanY) * 0.99, 0.01, -0.1,]); 
                    model = sprintf('a + b*exp(-1 * %010e * x)', c); % c specified up to 10 significant digits
                   ft = fittype(model, 'options', fo);
%                     ft = fittype('a + b*exp(c*x)', 'options', fo);
                    [fitobject, gof, output] = ...
                        fit(trialMeanX, trialMeanY, ft, fo);
%                     trialFit = fitobject.a + fitobject.b * exp(fitobject.c * x);      
                    trialFit = fitobject.a + fitobject.b * exp(-1 * c * x);      
%%

                    
                    Photometry.bleachFit(si, fCh).fitobject_trial = fitobject;
                    Photometry.bleachFit(si, fCh).gof_trial = gof;
                    Photometry.bleachFit(si, fCh).output_trial = output;
                    Photometry.bleachFit(si, fCh).trialTemplate = trialMeanY;                
                    Photometry.bleachFit(si, fCh).trialTemplateFull = trialMeanYFull;                    
                    Photometry.bleachFit(si, fCh).trialTemplateX = trialMeanX; % xData for session dFF fit
                    Photometry.bleachFit(si, fCh).trialTemplateFullX = (0:length(trialMeanYFull) - 1)' / sampleRate;
                    
                    Photometry.bleachFit(si, fCh).fitX =  x;
                    Photometry.bleachFit(si, fCh).trialFit = trialFit;                    
                    % set mean of true baseline period to zero
                    trialFit = trialFit - nanmean(trialFit(1, expFitStartP:blEndP));
                    blF = bsxfun(@plus, blF, trialFit);
                    dF = allData - blF;
                otherwise
            end            
            Photometry.data(fCh).dFF(tcounter:tcounter+nTrials - 1, :) = dF ./ blF; 
            Photometry.data(fCh).dF(tcounter:tcounter+nTrials - 1, :) = dF;
            sd = nanmean(nanstd(dF(:,blStartP:blEndP), 0, 2));
            Photometry.data(fCh).ZS(tcounter:tcounter+nTrials - 1, :) = dF / sd; % z-scored
            Photometry.data(fCh).raw(tcounter:tcounter+nTrials - 1, :) = allData;                  
            Photometry.data(fCh).blF(tcounter:tcounter+nTrials - 1, 1) = mean(blF, 2); 
            Photometry.data(fCh).blF_raw(tcounter:tcounter+nTrials - 1, 1) = blF_raw;            
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