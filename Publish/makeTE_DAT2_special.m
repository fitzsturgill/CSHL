function TE = makeTE_DAT2_special(sessions)

    if nargin < 1
        sessions = bpLoadSessions;
    end
    
    % find total number of trials acrosss selected sesssions
    scounter = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
    end
    nTrials = sum(scounter);
    statesToAdd= {'ITI', 'PreCsRecording', 'Cue', 'Delay', 'Reward', 'Punish', 'Omit', 'PostUsRecording'};

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
        'odorValve', zeros(nTrials, 1),...
        'ReinforcementOutcome', [],...
        'rewardSize', zeros(nTrials,1),...
        'rewardValveTime', zeros(nTrials, 1)...
        );
    
    TE.sessions = struct(... % sessions fields so that you can easily reload sessions data to modify TE
        'filename', cell(length(sessions), 1),...
        'filepath', [],...
        'index', [],...
        'OmitValveCode', []...
        );    

    for i = 1:length(statesToAdd)
        TE(1).(statesToAdd{i}) = bpAddStateAsTrialEvent(sessions, statesToAdd{i});
    end
    TE(1).Port1In = bpAddEventAsTrialEvent(sessions, 'Port1In');
    TE.filename = cell(nTrials, 1);
%     sessions(4).SessionData.TrialSettings(1).RewardValveTime
    tcounter = 1;
    for sCounter = 1:length(sessions)
        session = sessions(sCounter);
        TE.sessions(sCounter).filename = session.filename;
        TE.sessions(sCounter).filepath = session.filepath;
        TE.sessions(sCounter).index = sCounter;        
        if isfield(session.SessionData.Settings, 'NeturalToneOn')
            TE.sessions(sCounter).NeutralToneOn = session.SessionData.Settings.NeutralToneOn;
        else
            TE.sessions(sCounter).NeutralToneOn = 0;
        end
        
        if isfield(session.SessionData.Settings, 'OmitValveCode')
            TE.sessions(sCounter).OmitValveCode = session.SessionData.Settings.OmitValveCode;
        else
            TE.sessions(sCounter).OmitValveCode = 0;
        end        
            
        for counter = 1:session.SessionData.nTrials
            TE.filename{tcounter} = session.filename;
            TE.trialNumber(tcounter) = counter;
            TE.trialType(tcounter) = session.SessionData.TrialTypes(counter);
            TE.trialOutcome(tcounter) = session.SessionData.TrialOutcome(counter);
            TE.odorValve(tcounter) = session.SessionData.TrialSettings(counter).currentValve;
            if isfield(session.SessionData, 'Epoch')
                TE.epoch(tcounter) = session.SessionData.Epoch(counter);
            else
                TE.epoch(tcounter) = 1;
            end
            TE.sessionIndex(tcounter, 1) = sCounter;
            if any(isfinite(TE.Reward{tcounter}))
                TE.ReinforcementOutcome{tcounter} = 'Reward';
            elseif any(isfinite(TE.Punish{tcounter}))
                TE.ReinforcementOutcome{tcounter} = 'Punish';
            else
                TE.ReinforcementOutcome{tcounter} = 'Omit';                
            end
            TE.rewardValveTime(tcounter) = sessions(sCounter).SessionData.TrialSettings(counter).RewardValveTime;
            tcounter = tcounter + 1; % don't forget :)            
        end
    end
    
    % kludge for reward size
    smallTime = min(TE.rewardValveTime);
    bigTime = max(TE.rewardValveTime);
    TE.rewardSize(TE.rewardValveTime == smallTime) = 2;
    TE.rewardSize(TE.rewardValveTime == bigTime) = 8;
    TE.usZeros = cellfun(@(x) x(end), TE.Delay); %'Reward', 'Punish', 'WNoise', 'Neutral'
    
    TE.sessionChange = [0; diff(TE.sessionIndex)];    
    
    TE.cueCondition = ismember(TE.trialType, 1:2) * 1 + ismember(TE.trialType, 3:4) * 2; % 1 = rewarded odor, 2 = punished odor
    
   

    