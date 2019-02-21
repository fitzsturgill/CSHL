function Photometry = processTrialAnalysis_Photometry_expFitConcat(sessions, varargin)
% Fitz Sturgill 2018
% This function subtracts away a dual exponential bleaching trend by
% concatenating all trials together.  Assumes that bleaching occurs approximately without
% replenishment of unbleached fluorophore from a pool outside the
% photometry "field of view".  This assumtion is supported by the
% observation that fluorecense tends not to recover across successive days
% of photometry sessions.

% unlike other processTrialAnalysis_Photometry2, this function cuts away
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
    'plot', false;... %plots raw data and fit for quality control
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
maxSamplesPerSession = zeros(size(sessions));
for i = 1:length(sessions)
    scounter(i) = sessions(i).SessionData.nTrials;
    maxSamplesPerSession(i) = ceil(max(cellfun(@(x) size(x, 1), sessions(i).SessionData.NidaqData(:,1))) / s.downsample);
end
totalTrials = sum(scounter);

try
    sampleRate = sessions(1).SessionData.TrialSettings(1).nidaq.sample_rate;
catch
    sampleRate = 6100; % very early sessions don't have sample rate in settings
end
if s.uniformOutput % different sized trials are padded with NaNs, aligned to photometry start
    % determine maximum number of samples
    maxDuration = 0;
    for scounter = 1:length(sessions)
        for counter = 1:length(sessions(scounter).SessionData.TrialSettings)
            maxDuration = max(maxDuration, sessions(scounter).SessionData.TrialSettings(counter).nidaq.duration);
        end
    end
    originalSamples = maxDuration * sessions(1).SessionData.TrialSettings(1).nidaq.sample_rate;
    expectedSamples = ceil(originalSamples/s.downsample);
    if rem(sampleRate, s.downsample)
        error('downsample must be a factor of sampleRate');
    end
    % update sampleRate because you no longer need non-decimated sample rate
    sampleRate = sampleRate / s.downsample; % downsample must be a factor of sampleRate
    zeroTime = sessions(1).SessionData.RawEvents.Trial{1}.States.(s.zeroField)(1) - ...
        sessions(1).SessionData.RawEvents.Trial{1}.States.(s.startField)(1);
    xData = linspace(-zeroTime, ((expectedSamples - 1) / sampleRate) - zeroTime, expectedSamples);
else % unfinished, different sized trials are stored in cell arrays
    originalSamples = [];
    expectedSamples = [];
    sampleRate = sampleRate / s.downsample; % downsample must be a factor of sampleRate
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
        'blF_fit', NaN(totalTrials, expectedSamples),...
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
    if s.uniformOutput
        allData = NaN(nTrials, expectedSamples);
    else
        allData = NaN(nTrials, max(maxSamplesPerSession));
        allDataCell = cell(nTrials,1);
        sampleCount = zeros(nTrials, 1);
    end
    for chCounter = 1:length(s.channels)
        fCh = s.channels(chCounter);
        for trial = 1:nTrials
            trialData = SessionData.demod{trial, fCh}'; % convert to row vector
            if ismember(trial,sessions(si).excludeTrials)
                trialData(:)=nan;
            end

            if isempty(trialData)
                warning('Empty trial %d, think about it!',trial);
            end
            downData = decimate(trialData, s.downsample);
            downSamples = length(downData); % replace nSamples
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
                allDataCell{trial}= downData;
                sampleCount(trial) = length(downData);
            end
        end
        
        % baseline fluor.
        blStartP = bpX2pnt(s.baseline{chCounter}(1), sampleRate);
        blEndP = bpX2pnt(s.baseline{chCounter}(2), sampleRate);
        
        
        % concatenate all data together
        if s.uniformOutput
            blF_data = allData';
            blF_data = blF_data(:);
        else
            blF_data = cell2mat(allDataCell')';
        end
      	excludePoints=isnan(blF_data);
        fo = fitoptions('Method', 'NonlinearLeastSquares',...
            'Upper', [Inf range(blF_data(~excludePoints)) 0 range(blF_data(~excludePoints)) 0],...
            'Lower', [0 0 -2 0 -2],...    % 'Lower', [0 0 -1/5 0 -1/5],...
            'StartPoint', [min(blF_data(~excludePoints)) range(blF_data(~excludePoints))/2 -1 range(blF_data(~excludePoints))/2 -1]);
        model = 'a + b*exp(c*x) + d*exp(e*x)';
        ft = fittype(model, 'options', fo);
        x = (0:length(blF_data)-1)';
        [fitobject, gof, output] = ...
            fit(x(~excludePoints)-min(x(~excludePoints)), blF_data(~excludePoints), ft,fo);
        
        Photometry.bleachFit(si, fCh).fitobject_session = fitobject;
        Photometry.bleachFit(si, fCh).gof_session = gof;
        Photometry.bleachFit(si, fCh).output_session = output;
        blF_fit = nan(size(blF_data));
        blF_fit(~excludePoints) = feval(fitobject,x(~excludePoints)-min(x(~excludePoints)));
%         blF_fit=smoothdata(blF_data,'gaussian',1E4);
        
        if s.plot
            figure;
            plot(blF_data);
            hold on;
            %             blF_data_filter=(medfilt1(blF_data,3));
            %             plot([nan(3,1);blF_data_filter(4:end-4)]);
            blF_data_smooth=nan(size(blF_data));
            blF_data_smooth(~excludePoints)=smoothdata(blF_data(~excludePoints),'gaussian',1E4);
            plot(blF_data_smooth,'LineWidth',1.5);
            plot(blF_fit,'k','LineWidth',1.5);
            xlims=xlim;ylims=ylim;
            xtick=cumsum(sampleCount);
            xticklabel=1:length(sampleCount);
            xtickindex=1:1:length(xtick);
            set(gca,'XTick',xtick(xtickindex),'XTickLabel',xticklabel(xtickindex))
            text(xlims(1)+.1*diff(xlims),ylims(1)+.1*diff(ylims),sprintf('rsquare=%2.2f',gof.rsquare))
            title(sessions(si).filename,'Interpreter','none');
        end
        
        % subtract away the bleaching trend
        if s.uniformOutput
            blF_fit = reshape(blF_fit, size(allData, 2), size(allData, 1))';
            dF = allData - blF_fit;
        else
            blF_fit_reshaped=nan(size(allData));
            for t=1:length(sampleCount)
                blF_fit_reshaped(t,1:sampleCount(t))=blF_fit(1:sampleCount(t));
                blF_fit(1:sampleCount(t))=[];
            end
            blF_fit = blF_fit_reshaped;
            dF= allData - blF_fit;
        end
        if s.uniformOutput
            Photometry.data(fCh).dFF(tcounter:tcounter+nTrials - 1, :) = dF ./ blF_fit;
            Photometry.data(fCh).dF(tcounter:tcounter+nTrials - 1, :) = dF;
            sd = nanmean(nanstd(dF(:,blStartP:blEndP), 0, 2));
            Photometry.data(fCh).ZS(tcounter:tcounter+nTrials - 1, :) = dF / sd; % z-scored
            Photometry.data(fCh).raw(tcounter:tcounter+nTrials - 1, :) = allData;
            Photometry.data(fCh).blF_fit(tcounter:tcounter+nTrials - 1, :) = blF_fit;
            Photometry.data(fCh).ch = fCh;
        else
            dFF = dF ./ blF_fit;
            sd = nanmean(nanstd(dF(:,blStartP:blEndP), 0, 2));
            ZS = dF / sd;
            for trial = 1:nTrials
                Photometry.data(fCh).dFF{tcounter + trial - 1} = dFF(trial, 1:sampleCount(trial));
                Photometry.data(fCh).dF{tcounter + trial - 1} = dF(trial, 1:sampleCount(trial));
                Photometry.data(fCh).ZS{tcounter + trial - 1} = ZS(trial, 1:sampleCount(trial)); % z-scored
                Photometry.data(fCh).raw{tcounter + trial - 1} = allData(trial, 1:sampleCount(trial));
                Photometry.data(fCh).blF_fit{tcounter + trial - 1} = blF_fit(trial, 1:sampleCount(trial));
            end
    end

    end
    Photometry.startTime(tcounter:tcounter+nTrials - 1) = startTimes';
    tcounter = tcounter + nTrials;
    waitbar(si/length(sessions));
end
close(h);


