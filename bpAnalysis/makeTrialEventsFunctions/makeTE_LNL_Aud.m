function TE = makeTE_LNL_Aud(sessions)
    if nargin < 1
        sessions = bpLoadSessions;
    end
    
    % find total number of trials acrosss selected sesssions
    scounter = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
    end
    nTrials = sum(scounter);
    statesToAdd= {'NoLick', 'PreCsRecording', 'foreperiod', 'Cue', 'AnswerDelay', 'AnswerLick', 'AnswerNoLick',...
        'Reward', 'Punish', 'Neutral', 'PostUsRecording'};

    %% initialize TE
    TE = struct(...
        'filename', [],... 
        'trialNumber', zeros(nTrials, 1),...
        'trialType', zeros(nTrials, 1),...  
        'trialOutcome', NaN(nTrials, 1),... 
        'epoch', NaN(nTrials, 1),... 
        'sessionIndex', NaN(nTrials, 1),...
        'sessionChange', NaN(nTrials, 1),...
        'NoLick', zeros(nTrials, 1),...
        'Foreperiod', zeros(nTrials, 1),...
        'SoundValve', zeros(nTrials, 1),...
        'SoundValveIndex', zeros(nTrials, 1),...
        'CSValence', zeros(nTrials, 1),...
        'SoundAmplitude', zeros(nTrials, 1),...
        'BlockNumber', zeros(nTrials, 1),...                
        'BlockFcn', zeros(nTrials, 1),...                        
        'LickAction', [],...
        'Us', []...
        );

    for i = 1:length(statesToAdd)
        TE(1).(statesToAdd{i}) = bpAddStateAsTrialEvent(sessions, statesToAdd{i});
    end
    TE(1).Port1In = bpAddEventAsTrialEvent(sessions, 'Port1In');
    TE.filename = cell(nTrials, 1);
    TE.ReinforcementOutcome = cell(nTrials, 1);
    
        tcounter = 1;
    for sCounter = 1:length(sessions)
        session = sessions(sCounter);
        for counter = 1:session.SessionData.nTrials
            TE.filename{tcounter,1} = session.filename;
            TE.trialNumber(tcounter,1) = counter;
            TE.trialType(tcounter,1) = session.SessionData.TrialTypes(counter);
            TE.trialOutcome(tcounter,1) = session.SessionData.TrialOutcome(counter);
            TE.epoch(tcounter,1) = session.SessionData.Epoch(counter);
            TE.SoundValve(tcounter,1) = session.SessionData.SoundValve(counter);
            TE.SoundValveIndex(tcounter,1) = session.SessionData.SoundValveIndex(counter);
            TE.CSValence(tcounter,1) = session.SessionData.CSValence(counter);
            TE.BlockNumber(tcounter,1) = session.SessionData.BlockNumber(counter);
            TE.LickAction{tcounter,1} = session.SessionData.LickAction{counter};
            TE.ReinforcementOutcome{tcounter, 1} = session.SessionData.ReinforcementOutcome{counter};   
            TE.sessionIndex(tcounter, 1) = sCounter;
            usTimes = [TE.Reward{tcounter}; TE.Punish{tcounter}; TE.Neutral{tcounter}];
            TE.Us{tcounter, 1} = [max(usTimes(:,1)) max(usTimes(:,2))];
            TE.SoundAmplitude(tcounter,1)  = session.SessionData.SoundAmplitude(counter); 
            tcounter = tcounter + 1; % don't forget :)            
        end
    end
    

    TE.sessionChange = [0; diff(TE.sessionIndex)];
    TE.BlockChange = [0; diff(TE.BlockNumber) ~= 0];
    

    
%% wtf
nanOutcome = isnan(TE.trialOutcome);