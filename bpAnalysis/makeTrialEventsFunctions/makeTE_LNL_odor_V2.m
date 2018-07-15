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
        'filename', [],... 
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
        'BlockNumber', zeros(nTrials, 1),...                
        'BlockFcn', zeros(nTrials, 1),...                        
        'LickAction', [],...
        'Us', []...
        );
    
    % specific to auROC-driven reversals (uses auROC-driven block switch
    % function)
    if isfield(sessions(1).SessionData, 'AnswerLicksROC')
        rocOn = 1;
        TE.AnswerLicksROC = zeros(nTrials, 1);
    else
        rocOn = 0;
    end

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
            TE.OdorValve(tcounter,1) = session.SessionData.OdorValve(counter);
            TE.OdorValveIndex(tcounter,1) = session.SessionData.OdorValveIndex(counter);
            TE.CSValence(tcounter,1) = session.SessionData.CSValence(counter);
            TE.BlockNumber(tcounter,1) = session.SessionData.BlockNumber(counter);
            TE.LickAction{tcounter,1} = session.SessionData.LickAction{counter};
            TE.ReinforcementOutcome{tcounter, 1} = session.SessionData.ReinforcementOutcome{counter};   
            if rocOn
                TE.AnswerLicksROC(tcounter, 1) = session.SessionData.AnswerLicksROC.auROC(counter);
            end
            TE.sessionIndex(tcounter, 1) = sCounter;
            usTimes = [TE.Reward{tcounter}; TE.Punish{tcounter}; TE.WNoise{tcounter}; TE.Neutral{tcounter}];
            TE.Us{tcounter, 1} = [max(usTimes(:,1)) max(usTimes(:,2))];
            tcounter = tcounter + 1; % don't forget :)            
        end
    end
    

    TE.sessionChange = [0; diff(TE.sessionIndex)];
    TE.BlockChange = [0; diff(TE.BlockNumber) ~= 0];
    

    
%% wtf
nanOutcome = isnan(TE.trialOutcome);