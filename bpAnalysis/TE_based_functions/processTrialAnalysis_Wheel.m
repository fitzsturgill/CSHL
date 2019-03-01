function Wheel = processTrialAnalysis_Wheel(sessions, varargin)
% process wheel running speed measurements from wheel connected to rotary
% encoder.  Assumes that rotary pulses are unidirectional.
    %% optional parameters, first set defaults
    defaults = {...
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (store data in cell array)
        'startField', 'Start';...
        'zeroField', '';... % if empty, startField defines zero
        'Fs', 20;...
        'dpp', pi * 12.7 / 200; % distance per pulse
        'duration', 30;...
        'dataField', 'Port3In';... % field providing pulse times from rotary encoder
        'smoothWindow', 1;... % smoothing window in seconds
        };

    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings 
    if rem(s.duration * s.Fs, 1)
        error('*** duration must contain an integer number of samples ***');
    end    
    if isempty(s.zeroField)
        s.zeroField = s.startField;
    end
    % find total number of trials across selected sessions and size of
    % nidaq data
    scounter = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
    end
    totalTrials = sum(scounter);        
    if s.uniformOutput        
        zeroTime = sessions(1).SessionData.RawEvents.Trial{1}.States.(s.zeroField)(1) - ...
            sessions(1).SessionData.RawEvents.Trial{1}.States.(s.startField)(1);
        xData = -zeroTime:1/s.Fs:(-zeroTime + s.duration) - 1/s.Fs;
    else
        error(' uniformOutput == 0 mode not yet implemented ');
    end
    
    Wheel = struct(...
        'data', [],...
        'settings', s,...
        'startTime', NaN(totalTrials, 1),...
        'Fs', s.Fs,...
        'xData', xData...
        );
    s.smoothWindow = s.smoothWindow * s.Fs; % convert smooth window from time to samples
    if s.uniformOutput
        data = struct(...
            'V', NaN(totalTrials, length(xData)),...
            'X', NaN(totalTrials, length(xData))...            
            );
    else
        data = struct(...
            'V', {},...
            'X', {}...
            );        
    end
    

    Wheel.data = data;
    tcounter = 1;    
    for si = 1:length(sessions)
        SessionData = sessions(si).SessionData;
        nTrials = SessionData.nTrials;
        for trial = 1:nTrials
            startTime = SessionData.RawEvents.Trial{trial}.States.(s.startField)(1);
            if isfield(SessionData.RawEvents.Trial{trial}.Events, s.dataField)
                pulseTimes = SessionData.RawEvents.Trial{trial}.Events.(s.dataField) - startTime;
            else
                pulseTimes = [];
            end
            edges = 0:1/s.Fs:s.duration;
            position = cumsum(histcounts(pulseTimes, edges));
            try
                position = smooth(position, s.smoothWindow);
            catch
                position = smooth(position, 'linear', s.smoothWindow);
            end
            velocity = gradient(position); % gradient preserves number of points unlike diff
            Wheel.startTime(tcounter) = startTime;
            Wheel.data.X(tcounter,:) = position;
            Wheel.data.V(tcounter,:) = velocity;            
            tcounter = tcounter + 1;              
        end                                          
    end
                