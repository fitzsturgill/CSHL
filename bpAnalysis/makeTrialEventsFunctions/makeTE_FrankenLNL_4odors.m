function TE = makeTE_FrankenLNL_4odors(sessions)
    if nargin < 1
        sessions = bpLoadSessions;
    end
    
    % find total number of trials acrosss selected sesssions
    sTally = zeros(size(sessions));
    for i = 1:length(sessions)
        sTally(i) = sessions(i).SessionData.nTrials;
    end
    nTrials = sum(sTally);
    statesToAdd= {'ITI', 'PreCsRecording', 'Cue1', 'Trace1', 'Cue2', 'Trace2', 'Outcome', 'Reward', 'Punish', 'WNoise', 'Neutral', 'PostUsRecording'};

    %% initialize TE
    TE = struct(...
        'filename', [],... 
        'trialNumber', zeros(nTrials, 1),...
        'trialType', zeros(nTrials, 1),...  
        'trialOutcome', NaN(nTrials, 1),... 
        'sessionIndex', NaN(nTrials, 1),...
        'sessionChange', NaN(nTrials, 1),...
        'ITI', zeros(nTrials, 1),...
        'Odor1Valve', zeros(nTrials, 1),...
        'Odor1ValveIndex', zeros(nTrials, 1),...
        'Odor2Valve', zeros(nTrials, 1),...
        'Odor2ValveIndex', zeros(nTrials, 1),...        
        'BlockNumber', zeros(nTrials, 1),...                
        'BlockFcn', zeros(nTrials, 1),...                        
        'Us', [],...
        'TrialStartTimestamp', [],...
        'sessions', []... % new 9/2018, store names and filepaths and sessionIndex, this struct is of length nsessions with fields: names, filespaths, indices
        );
    
    TE.sessions = struct(... % sessions fields so that you can easily reload sessions data to modify TE
        'filename', cell(length(sessions), 1),...
        'filepath', [],...
        'index', []...
        );
    
    %% updated for new version of FrankenLNL with tone and light options in blocks           
    additionalFields = {'CS1_tone', 'CS2_tone', 'CS1_light', 'CS2_light'}; 
    counter = 1;
    while length(additionalFields) >= counter
        if isfield(sessions(1).SessionData, additionalFields{counter})
            TE.(additionalFields{counter}) = zeros(nTrials, 1);
            counter = counter + 1;
        else
            keep = true(length(additionalFields), 1);
            keep(counter) = false;
            additionalFields = additionalFields(keep);
        end
    end
    %%
    for i = 1:length(statesToAdd)
        TE(1).(statesToAdd{i}) = bpAddStateAsTrialEvent(sessions, statesToAdd{i});
    end
    TE(1).Port1In = bpAddEventAsTrialEvent(sessions, 'Port1In');
    TE.filename = cell(nTrials, 1);
    TE.ReinforcementOutcome = cell(nTrials, 1);
    
        tcounter = 1;
    for sCounter = 1:length(sessions)
        session = sessions(sCounter);
        TE.sessions(sCounter).filename = session.filename;
        TE.sessions(sCounter).filepath = session.filepath;
        TE.sessions(sCounter).index = sCounter;
        for counter = 1:session.SessionData.nTrials
            TE.filename{tcounter,1} = session.filename;
            TE.trialNumber(tcounter,1) = counter;
            TE.trialType(tcounter,1) = session.SessionData.TrialTypes(counter);
            TE.Odor1Valve(tcounter,1) = session.SessionData.Odor1Valve(counter);
            TE.Odor1ValveIndex(tcounter,1) = session.SessionData.Odor1ValveIndex(counter);
            TE.Odor2Valve(tcounter,1) = session.SessionData.Odor2Valve(counter);
            TE.Odor2ValveIndex(tcounter,1) = session.SessionData.Odor2ValveIndex(counter);            
            TE.BlockNumber(tcounter,1) = session.SessionData.BlockNumber(counter);
            TE.ReinforcementOutcome{tcounter, 1} = session.SessionData.ReinforcementOutcome{counter};
            TE.TrialStartTimestamp(tcounter, 1) = session.SessionData.TrialStartTimestamp(counter);
            TE.sessionIndex(tcounter, 1) = sCounter;
            usTimes = [TE.Reward{tcounter}; TE.Punish{tcounter}; TE.WNoise{tcounter}; TE.Neutral{tcounter}];
            TE.Us{tcounter, 1} = [max(usTimes(:,1)) max(usTimes(:,2))];
            for afCounter = 1:length(additionalFields)
                TE.(additionalFields{afCounter})(tcounter, 1) = session.SessionData.(additionalFields{afCounter})(counter);
            end
            tcounter = tcounter + 1; % don't forget :)            
        end
    end
    

    TE.sessionChange = [0; diff(TE.sessionIndex)];
    TE.BlockChange = [0; diff(TE.BlockNumber) ~= 0];
      