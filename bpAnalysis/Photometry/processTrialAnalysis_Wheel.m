function Photometry = processTrialAnalysis_Wheel(sessions, varargin)

% Fs - sample rate
% Ns - number of samples/trial  (currently a scalar value, uniformOutput =
% 0 is not yet implemented
    
% !!!! Note you can have wheel increment values outside of the time span of
% the photometry acquisition, you need startField or zero or something, dT and nSamples
    %% optional parameters, first set defaults
    defaults = {...
        'dX', 5;... %unit length in cm traveled for each advance of the rotary encoder, (currently a weakly founded guess)
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (store data in cell array)
        'startField', 'Start';...
        'zeroField', '';... % if empty, startField defines zero
        'Fs', 20;...
        'duration', 30;...
        'dataField', 'Port3In';...
        };
    
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings 
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
        xData = linspace(-zeroTime, -zeroTime + s.duration, round(s.duration/s.Fs));
    else
        error(' uniformOutput == 0 mode not yet implemented ');
    end
    
    Wheel = struct(...
        'data', [],...
        'settings', s,...
        'startTime', NaN(totalTrials, 1),...
        'xData', xData...
        );
    if s.uniformOutput
        data = struct(...
            'V', NaN(totalTrials, length(xData)),...
            );
    else
        data = struct(...
            'V', {}...
            );        
    end
    
     %%
D = 12.7; % diameter of wheel in cm
pr = 200; % pulses/rotation 
dpp = pi * 12.7 / 200;
rawTimes = sessions(1).SessionData.RawEvents.Trial{5}.Events.Port3In;
wheelTimes = rawTimes(rawTimes < 30);
edges = 0:1/20:30;
% wheel_X = cumsum(histcounts(rawTimes, edges) * dpp);
wheel_X = smooth(cumsum(histcounts(rawTimes, edges)), 20);  
wheel_t = 0:1/20:30-1/20;
  
ensureFigure('test', 1);
subplot(2,1,1);
plot(wheelTimes, cumsum(ones(size(wheelTimes)) * dpp), 'b'); hold on;
plot(wheel_t, wheel_X, 'g');

subplot(2,1,2);
plot(wheel_t(1:end-1), diff(wheel_X), 'b'); hold on
plot(wheel_t, gradient(wheel_X), 'g');
 
    tcounter = 1;    
    for si = 1:length(sessions);
        SessionData = sessions(si).SessionData;
        startTimes = cellfun(@(x) x.States.(s.startField)(1), sessions(si).SessionData.RawEvents.Trial); % take the beginning time stamp for the startField-specified Bpod state
        nTrials = SessionData.nTrials;
        allData = NaN(nTrials, newSamples);   
            for trial = 1:nTrials
                
                
            end
        Wheel.startTime(tcounter:tcounter+nTrials - 1) = startTimes';
        tcounter = tcounter + nTrials;            
    end
                