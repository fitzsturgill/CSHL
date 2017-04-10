function TE = makeTE_LNL_odor_V2(sessions)
    if nargin < 1
        sessions = bpLoadSessions;
    end
    
    % find total number of trials acrosss selected sesssions
    scounter = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
    end
    nTrials = sum(scounter);
    statesToAdd= {'ITI', 'PreCsRecording', 'Cue', 'AnswerDelay', 'AnswerStart', 'AnswerLick', 'AnswerNoLick',...
        'Reward', 'Punish', 'WNoise', 'Neutral', 'PostUsRecording'};

    %% initialize TE
    TE = struct(...
        'filename', {},... 
        'trialNumber', zeros(nTrials, 1),...
        'trialType', zeros(nTrials, 1),...  
        'trialOutcome', NaN(nTrials, 1),... 
        'epoch', NaN(nTrials, 1),... 
        'sessionIndex', NaN(nTrials, 1),...
        'sessionChange', NaN(nTrials, 1),...
        'ITI', zeros(nTrials, 1),...
        'OdorValve', zeros(nTrials, 1),...
        'OdorValveIndex', zeros(nTrials, 1),...
        'CSValence', zeros(nTrials, 1),...        
        'ReinforcementOutcome', {},...            
        'BlockNumber', zeros(nTrials, 1),...                
        'BlockFcn', zeros(nTrials, 1),...                        
        'LickAction', {}...                        
        );

    for i = 1:length(statesToAdd)
        TE(1).(statesToAdd{i}) = bpAddStateAsTrialEvent(sessions, statesToAdd{i});
    end
    TE(1).Port1In = bpAddEventAsTrialEvent(sessions, 'Port1In');
    TE.filename = cell(nTrials, 1);
    
        tcounter = 1;
    for sCounter = 1:length(sessions)
        session = sessions(sCounter);
        for counter = 1:session.SessionData.nTrials
            TE.filename{tcounter} = session.filename;
            TE.trialNumber(tcounter) = counter;
            TE.trialType(tcounter) = session.SessionData.TrialTypes(counter);
            TE.trialOutcome(tcounter) = session.SessionData.TrialOutcome(counter);
            TE.epoch(tcounter) = session.SessionData.Epoch(counter);
            TE.OdorValve(tcounter) = session.SessionData.OdorValve(counter);
            TE.OdorValveIndex(tcounter) = session.SessionData.OdorValveIndex(counter);
            TE.CSValence(tcounter) = session.SessionData.CSValence(counter);
            TE.BlockNumber(tcounter) = session.SessionData.BlockNumber(counter);
            TE.LickAction{tcounter} = session.SessionData.LickAction{counter};
            tcounter = tcounter + 1; % don't forget :)            
        end
    end
    
    sessionNames = unique(TE.filename);
    for counter = 1:length(sessionNames)
        sname = sessionNames{counter};
        TE.sessionIndex(cellfun(@(x) strcmp(x, sname), TE.filename)) = counter;
    end
    TE.sessionIndex = TE.sessionIndex';
    TE.sessionChange = [0; diff(TE.sessionIndex)];
    
    
    
