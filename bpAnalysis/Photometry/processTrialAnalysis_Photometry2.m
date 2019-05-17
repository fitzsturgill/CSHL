function Photometry = processTrialAnalysis_Photometry2(sessions, varargin)
% exemplar for new trial analysis functions


    %% optional parameters, first set defaults
    defaults = {...
        'channels', 1;... % 8/28/2016- changed channels default from [] to 1
        'refChannels', [];...
        'baseline', [1 4];... % 1 - 3 second into recordinge
        'blMode', 'byTrial';... % 'byTrial', 'bySession', 'expFit'  expFit- interpolates baseline from biexponential fit to raw fl baselines across trials
        'dFFMode', 'simple';... % 'simple', 'expFit'   !(now with hard-coded time constant for exponential) expFit- subtracts within-trial exponential bleaching trend using an exponential fit to the trial average baseline period
        'expFitBegin', 0.1;...
        'zeroField', 'Us';...
        'startField', 'PreCsRecording';... % TO DO: PROVIDE AUTOMATICALLY BY BPOD NIDAQ CODE 
        'downsample', 305;...  % much faster, relies upon decimation, originally specified smaple rate must be evenly divisible by decimation factor
        'sampleRate', []; % often much slower, use to achieve an exact final sample rate
        'uniformOutput', 1;...            % set to 0 if acqs are variable in length (store data in cell array)
%         'tau', [];... % tau option currently deprecated, now tau is a
%         free parameter
        'forceAmp', 0;... % % force demodulation even if the refChannel LED is off (i.e. it's amplitude = 0)
        'ACfilter', 0;...
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    %% if not specified per channel (as a cell array), use identical values across channels
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
    
    
    % get old sample rate    
    try
        if isempty(s.sampleRate) % for downsample/decimation, use requested not achieved nidaq sample rate
            sampleRate_old = sessions(1).SessionData.TrialSettings(1).nidaq.sample_rate;
        else % to achieve an exact resampled sample rate, try to get the achieved nidaq sample rate (a la setverify)
            try
                sampleRate_old = sessions(1).SessionData.NidaqData{1,2}.sample_rate(1); 
            catch
                sampleRate_old = sessions(1).SessionData.TrialSettings(1).nidaq.sample_rate;    
            end
        end
    catch
        sampleRate_old = 6100; % very early sessions don't have sample rate in settings
    end
    
    
    if ~isempty(s.sampleRate)% final sample rate and coefficients for resampling    
        sampleRate = s.sampleRate;
        [p, q] = rat(sampleRate/sampleRate_old);            
    else 
        if rem(sampleRate_old, s.downsample)
            error('downsample must be a factor of sampleRate');
        end
        sampleRate = sampleRate_old / s.downsample;
    end

    
    
    % find total number of trials across selected sessions and maximum size of
    % nidaq data per session
    scounter = zeros(size(sessions));    
    maxSamplesPerSession = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
        if ~isempty(s.sampleRate)
            maxSamplesPerSession(i) = ceil(max(cellfun(@(x) ceil(size(x, 1) * p/q), sessions(i).SessionData.NidaqData(:,1))));
        else
            maxSamplesPerSession(i) = ceil(max(cellfun(@(x) ceil(size(x, 1) / s.downsample), sessions(i).SessionData.NidaqData(:,1))));
        end
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
        originalSamples = ceil(maxDuration * sampleRate_old);
        if ~isempty(s.sampleRate)
            expectedSamples = ceil(originalSamples * p/q);                   
        else
            expectedSamples = ceil(originalSamples / s.downsample);  
        end
        zeroTime = sessions(1).SessionData.RawEvents.Trial{1}.States.(s.zeroField)(1) - ...
            sessions(1).SessionData.RawEvents.Trial{1}.States.(s.startField)(1);
        xData = linspace(-zeroTime, ((expectedSamples - 1) / sampleRate) - zeroTime, expectedSamples);    
    else % different sized trials are stored in cell arrays, no padding occurs
        originalSamples = [];
        expectedSamples = [];
        xData = [];
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
            'dFF', NaN(totalTrials, expectedSamples),... % deltaF/F
            'dF', NaN(totalTrials, expectedSamples),... % deltaF            
            'ZS', NaN(totalTrials, expectedSamples),... % deltaF/F            
            'raw', NaN(totalTrials, expectedSamples),... %
            'blF', NaN(totalTrials, 1),...
            'blF_raw', NaN(totalTrials, 1),...
            'blF_fit', NaN(totalTrials, 1),...
            'ch', []...
            );    
    else % unfinished, different sized trials are stored in cell arrays
        data = struct(...        
            'dFF', {},... % deltaF/F
            'dF', {},...
            'ZS', {},... % deltaF/F            
            'raw', {},... %
            'blF', NaN(totalTrials, 1),...
            'blF_raw', NaN(totalTrials, 1),...
            'blF_fit', NaN(totalTrials, 1),...
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

    
    tcounter = 1;    
    for si = 1:length(sessions)
        SessionData = sessions(si).SessionData;
        if ~isfield(SessionData, 'demod')
            SessionData = demodulateSession(SessionData, 'channels', s.channels, 'refChannels', s.refChannels, 'forceAmp', s.forceAmp, 'ACfilter', s.ACfilter); % don't necessarily want to save these back to sessions because that'd eat up memory
        end
        startTimes = cellfun(@(x) x.States.(s.startField)(1), sessions(si).SessionData.RawEvents.Trial); % take the beginning time stamp for the startField-specified Bpod state
        nTrials = SessionData.nTrials;
        if s.uniformOutput
            allData = NaN(nTrials, expectedSamples);   
        else
            allData = NaN(nTrials, maxSamplesPerSession(si));
            sampleCount = zeros(nTrials, 1);
        end


        for chCounter = 1:length(s.channels)
            fCh = s.channels(chCounter);
            blMode = s.blMode{chCounter};
            dFFMode = s.dFFMode{chCounter};
%             blF_fit = NaN(nTrials, 1);
%             blF_raw = NaN(nTrials, 1);


            
            for trial = 1:nTrials
                trialData = SessionData.demod{trial, fCh}'; % convert to row vector

                if isempty(trialData)
                    warning('do something for empty trials here, maybe continue?');
                end
               
                if isempty(s.downsample)
                    downData = resample(trialData, p, q);
                else
                    downData = decimate(trialData, s.downsample);
                end
                downSamples = length(downData); % replace nSamples
                
                % baseline fluor.
                blStartP = min(bpX2pnt(s.baseline{chCounter}(1), sampleRate), downSamples);
                blEndP = min(bpX2pnt(s.baseline{chCounter}(2), sampleRate), downSamples); 
                expFitStartP = min(bpX2pnt(s.expFitBegin{chCounter}, sampleRate), downSamples);
                
                if s.uniformOutput
                % in case nidaq acquisition ended early for some reason, pad with NaNs, this should be fixed as of 8/2016                
                    if downSamples < expectedSamples
                        downData = [downData NaN(1, expectedSamples - downSamples)];
                    elseif downSamples > expectedSamples
                        downData = downData(1:expectedSamples);
                    end               
                    allData(trial, :) = downData;
                else
                    allData(trial, 1:downSamples) = downData;
                    sampleCount(trial) = length(downData);
                end
            end


            
            
            % always perform session bleaching fit and store it for later
            % use

            blF_fit = nanmean(allData(:, expFitStartP:blEndP), 2); % take mean across time, not trials                   

            fo = fitoptions('Method', 'NonlinearLeastSquares',...
                'Upper', [Inf range(blF_fit) 0 range(blF_fit) 0],...
                'Lower', [0 0 -1 0 -1],...    % 'Lower', [0 0 -1/5 0 -1/5],...                    
                'StartPoint', [min(blF_fit) range(blF_fit)/2 -5 range(blF_fit)/2 -100]...
                );
            model = 'a + b*exp(c*x) + d*exp(e*x)';
            ft = fittype(model, 'options', fo);
            if any(isnan(blF_fit))
                blF_fit = inpaint_nans(blF_fit);
            end
            x = (0:length(blF_fit)-1)';
%             [fitobject, gof, output] = ...  % old
%                 fit((1:size(allData, 1))', blF_fit, ft, fo);
            [fitobject, gof, output] = ... % new
                fit(x, blF_fit, ft, fo);
%             warning('check this in processTrialAnalysis_Photometry2, x should start at 0 for session bleach fit'); % old --> new
            Photometry.bleachFit(si, fCh).fitobject_session = fitobject;
            Photometry.bleachFit(si, fCh).gof_session = gof;
            Photometry.bleachFit(si, fCh).output_session = output;

            blF_fit = fitobject.a + fitobject.b * exp(fitobject.c * x) + fitobject.d * exp(fitobject.e * x);   
            
            % apply the selected baseline mode

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
                    blF = repmat(blF_fit, 1, size(allData, 2));
                    blF_raw = nanmean(allData(:, blStartP:blEndP), 2);
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

%% single exponential, 
                    fo = fitoptions('Method', 'NonlinearLeastSquares',...
                        'Lower', [0 0 0.5],...%, -Inf],...
                        'Upper', [max(trialMeanY) range(trialMeanY) 4],... % old upper tau = 4
                        'StartPoint', [min(trialMeanY) * 0.99, 0.01 2]);%, -0.3,]);

                    model = 'a + b*exp(-1 /c * x)';
                    ft = fittype(model, 'options', fo);
                    [fitobject, gof, output] = ...
                        fit(trialMeanX, trialMeanY, ft, fo);
                    trialFit = fitobject.a + fitobject.b * exp(-1 /fitobject.c * x);      
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
            if s.uniformOutput
                Photometry.data(fCh).dFF(tcounter:tcounter+nTrials - 1, :) = dF ./ blF; 
                Photometry.data(fCh).dF(tcounter:tcounter+nTrials - 1, :) = dF;
                sd = nanmean(nanstd(dF(:,blStartP:blEndP), 0, 2));
                Photometry.data(fCh).ZS(tcounter:tcounter+nTrials - 1, :) = dF / sd; % z-scored
                Photometry.data(fCh).raw(tcounter:tcounter+nTrials - 1, :) = allData;                  
                Photometry.data(fCh).blF(tcounter:tcounter+nTrials - 1, 1) = mean(blF, 2); 

            else
                dFF = dF ./ blF;
                sd = nanmean(nanstd(dF(:,blStartP:blEndP), 0, 2));
                ZS = dF / sd;
                blF_vector = nanmean(blF, 2);
                for trial = 1:nTrials
                    Photometry.data(fCh).dFF{tcounter + trial - 1} = dFF(trial, 1:sampleCount(trial)); 
                    Photometry.data(fCh).dF{tcounter + trial - 1} = dF(trial, 1:sampleCount(trial));
                    Photometry.data(fCh).ZS{tcounter + trial - 1} = ZS(trial, 1:sampleCount(trial)); % z-scored
                    Photometry.data(fCh).raw{tcounter + trial - 1} = allData(trial, 1:sampleCount(trial));                  
                end
            end
            Photometry.data(fCh).blF(tcounter:tcounter+nTrials - 1, 1) = nanmean(blF, 2);                 
            Photometry.data(fCh).blF_raw(tcounter:tcounter+nTrials - 1, 1) = blF_raw;         
            Photometry.data(fCh).blF_fit(tcounter:tcounter+nTrials - 1, 1) = blF_fit;
            Photometry.data(fCh).ch = fCh;
        end
        Photometry.startTime(tcounter:tcounter+nTrials - 1) = startTimes';
        tcounter = tcounter + nTrials;
        waitbar(si/length(sessions));
    end
    close(h);
    

    
%         originalSamples = zeros(1,length(sessions));
%         for si = 1:length(sessions);
%             originalSamples(si) = max(cellfun(@(x) size(x,1), sessions(si).SessionData.NidaqData(:,1)));
%         end
%         originalSamples = max(originalSamples);
