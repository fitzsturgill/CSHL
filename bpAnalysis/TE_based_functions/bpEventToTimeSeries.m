function eventTS = bpEventToTimeSeries(TE, event, varargin)

    %% optional parameters, first set defaults
    defaults = {...
        'zeroField', 'Us';...
        'startField', 'PreCsRecording';... %
        'sampleRate', 20;...
        'uniformOutput', 1;...            % not currently implemented, idea is to set to 0 if acqs are going to be variable in length (store data in cell array)
        'duration', 11;...
        };
    [s, ~] = parse_args(defaults, varargin{:}); % combine default and passed (via varargin) parameter settings
    
    nTrials = length(TE.filename);

    %% following code is kind of weak, just following along with functions like processTrialAnalysis_Photometry2 which are different in that they take as their input sessions structures
    if s.uniformOutput        
        zeroTime = TE.(s.zeroField){1}(1) - ... % assume that first trial is indicative of rest, i.e. uniformOutput  = true
            TE.(s.startField){1}(1);
        xData = -zeroTime:1/s.sampleRate:(-zeroTime + s.duration) - 1/s.sampleRate;
    else
        error(' uniformOutput == 0 mode not yet implemented ');
    end
    
    %%
    eventTS = struct(...
        'settings', s,...
        'data', zeros(nTrials, length(xData)),...
        'xData', xData...
        );

    
    rawEdges = 0:1/s.sampleRate:s.duration; % length should be 1 greater than xData since result of histcounts is of length one less than edges
    for trial = 1:nTrials
        startTime = TE.(s.startField){trial}(1);
        edges = startTime + rawEdges;
        
        theseEvents = TE.(event){trial};
        if ~isnan(theseEvents(1))
            eventTS.data(trial, :) = histcounts(theseEvents, edges);
        end
    end