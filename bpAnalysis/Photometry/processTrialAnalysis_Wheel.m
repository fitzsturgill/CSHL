function Photometry = processTrialAnalysis_Wheel(sessions, Fs, Ns, varargin)

% Fs - sample rate
% Ns - number of samples/trial  (currently a scalar value, uniformOutput =
% 0 is not yet implemented
    
% !!!! Note you can have wheel increment values outside of the time span of
% the photometry acquisition, you need startField or zero or something, dT and nSamples
    %% optional parameters, first set defaults
    defaults = {...
        'dX', 5;... %unit length in cm traveled for each advance of the rotary encoder, (currently a weakly founded guess)
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (store data in cell array)
        'startField', 'PreCsRecording'
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    s.Fs = Fs;    
    
    % find total number of trials across selected sessions and size of
    % nidaq data
    scounter = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
    end
    totalTrials = sum(scounter);    

    if s.uniformOutput 

        xData = linspace(-zeroTime, ((newSamples - 1) / sampleRate) - zeroTime, newSamples);    