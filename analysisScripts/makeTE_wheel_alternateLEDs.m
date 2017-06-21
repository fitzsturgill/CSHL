function TE = makeTE_wheel_alternateLEDs(sessions)
    if nargin < 1
        sessions = bpLoadSessions;
    end
    
    % find total number of trials acrosss selected sesssions
    scounter = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
    end
    nTrials = sum(scounter);

    %% initialize TE
    TE = struct(...
        'filename', [],... 
        'trialNumber', zeros(nTrials, 1),...
        'trialType', zeros(nTrials, 1),...  
        'trialStartTimeStamp', zeros(nTrials, 1),...
        'Epoch', NaN(nTrials, 1),... 
        'sessionIndex', NaN(nTrials, 1),...
        'sessionChange', NaN(nTrials, 1),...
        'ITI', zeros(nTrials, 1),...
        'LED1_amp', zeros(nTrials, 1),...
        'LED2_amp', zeros(nTrials, 1)...        
        );

    TE(1).Baseline = bpAddStateAsTrialEvent(sessions, 'Baseline');
    TE(1).Reward = bpAddStateAsTrialEvent(sessions, 'Reward*');
    TE(1).Port1In = bpAddEventAsTrialEvent(sessions, 'Port1In');
    TE.filename = cell(nTrials, 1);
    
    tcounter = 1;
    for sCounter = 1:length(sessions)
        session = sessions(sCounter);
        for counter = 1:session.SessionData.nTrials
            TE.filename{tcounter,1} = session.filename;
            TE.trialNumber(tcounter,1) = counter;
            TE.trialStartTimeStamp(tcounter) = session.SessionData.TrialStartTimestamp(counter);
            if isfield(session.SessionData, 'Epoch')
                TE.Epoch(tcounter,1) = session.SessionData.Epoch(counter);
            end
            TE.LED1_amp(tcounter) = session.SessionData.NidaqData{counter, 2}.amp(1);
            TE.LED2_amp(tcounter) = session.SessionData.NidaqData{counter, 2}.amp(2);            
            tcounter = tcounter + 1; % don't forget :)    
        end
    end
    
    sessionNames = unique(TE.filename);
    for counter = 1:length(sessionNames)
        sname = sessionNames{counter};
        TE.sessionIndex(cellfun(@(x) strcmp(x, sname), TE.filename)) = counter;
    end
    TE.sessionChange = [0; diff(TE.sessionIndex)];
    
    
    
