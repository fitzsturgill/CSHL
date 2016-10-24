function TE = makeTE_CuedOutcome_Odor_Complete(sessions)
    if nargin < 1
        sessions = bpLoadSessions;
    end
    
    % find total number of trials acrosss selected sesssions
    scounter = zeros(size(sessions));
    for i = 1:length(sessions)
        scounter(i) = sessions(i).SessionData.nTrials;
    end
    nTrials = sum(scounter);
    statesToAdd= {'ITI', 'PreCsRecording', 'Cue', 'Delay', 'Us', 'PostUsRecording'};

    %% initialize TE
    TE = struct(...
        'filename', {},... 
        'trialNumber', zeros(nTrials, 1),...
        'trialType', zeros(nTrials, 1),...  
        'trialOutcome', NaN(nTrials, 1),... 
        'epoch', NaN(nTrials, 1),... 
        'sessionIndex', NaN(nTrials, 1),...
        'sessionChange', NaN(nTrials, 1),...
        'ITI', zeros(nTrials, 1)...
        );

    for i = 1:length(statesToAdd)
        TE(1).(statesToAdd{i}) = bpAddStateAsTrialEvent(sessions, statesToAdd{i});
    end
    TE(1).Port1In = bpAddEventAsTrialEvent(sessions, 'Port1In');
    TE.filename = cell(nTrials, 1);

    %% bpCountEventByStates2 is outdated
%     TE.csLicks = bpCountEventByStates2(sessions, 'Port1In', 'Cue', 'endState', 'Delay'); % #licks spanning Cue and Delay
%     TE.usLicks = bpCountEventByStates2(sessions, 'Port1In', 'Us', 'endState', 'PostUsRecording'); % #licks spanning Us and PostUsRecording

    
    tcounter = 1;
    for sCounter = 1:length(sessions)
        session = sessions(sCounter);
        for counter = 1:session.SessionData.nTrials
            TE.filename{tcounter} = session.filename;
            TE.trialNumber(tcounter) = counter;
            TE.trialType(tcounter) = session.SessionData.TrialTypes(counter);
            TE.trialOutcome(tcounter) = session.SessionData.TrialOutcome(counter);
            TE.epoch(tcounter) = session.SessionData.Epoch(counter);
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
    
    
    
%     
%     function TE = makeTE_CuedOutcome_Odor_Complete(TE, sessions)
%     if nargin < 2    
%         sessions = bpLoadSessions;
%     end
%     %% make structure with many different fields, e.g. trialNumber, trialOutcome, reactiontime, etc.
%     % every field in this structure should be either a numeric array with numRows
%     % (dimension 1) == nTrials OR a cell array with the same criteria
% 
%     % find total number of trials acrosss selected sesssions
%     scounter = zeros(size(sessions));
%     for i = 1:length(sessions)
%         scounter(i) = sessions(i).SessionData.nTrials;
%     end
%     
%     nTrials = sum(scounter);
%     statesToAdd= {'ITI', 'PreCsRecording', 'Cue', 'Delay', 'Us', 'PostUsRecording'};
%     if nargin < 1 || isempty(TE) % initialize from scratch (reload)
%         startIndex = 1;
%         %by default, TE... filled with NaN unless conditions are met
%         TE = struct(...
%             'filename', {},... 
%             'trialNumber', zeros(nTrials, 1),...
%             'trialType', zeros(nTrials, 1),...  
%             'trialOutcome', NaN(nTrials, 1),... 
%             'epoch', NaN(nTrials, 1),... 
%             'csLicks', {},...
%             'usLicks', {},...
%             'sessionIndex', NaN(nTrials, 1),...
%             'sessionChange', NaN(nTrials, 1),...
%             'ITI', zeros(nTrials, 1)...
%             );
% 
%             for i = 1:length(statesToAdd)
%                 TE(1).(statesToAdd{i}) = bpAddStateAsTrialEvent(sessions, statesToAdd{i});
%             end
%             TE(1).Port1In = bpAddEventAsTrialEvent(sessions, 'Port1In');
%             TE.filename = cell(nTrials, 1);
%     else
%         startIndex = length(TE.filename);
%         % extend all fields from startIndex
%         TE.filename{startIndex + nTrials} = ''; 
%         TE.trialNumber = [TE.trialNumber ; zeros(nTrials, 1)];
%         TE.trialType = [TE.trialType; zeros(nTrials, 1)];  
%         TE.trialOutcome = [TE.trialOutcome; NaN(nTrials, 1)];  
%         TE.epoch = [TE.epoch; zeros(nTrials, 1)]; 
%         TE.csLicks{startIndex + nTrials} = [];
%         TE.usLicks{startIndex + nTrials} = [];        
%         TE.sessionIndex = [TE.sessionIndex; NaN(nTrials, 1)];
%         TE.sessionChange = [TE.sessionChange; NaN(nTrials, 1)];
%         TE.ITI= [TE.ITI; zeros(nTrials, 1)];
%         for i = 1:length(statesToAdd)
%             TE.(statesToAdd{i}) = bpAddStateAsTrialEvent(sessions, statesToAdd{i}, TE.(statesToAdd{i}));
%         end
%         TE.Port1In = bpAddEventAsTrialEvent(sessions, 'Port1In', TE.Port1In);
%         TE.filename = cell(nTrials, 1);        
%     end
% 
%     TE.csLicks = 
%     tcounter = 1;
%     for sCounter = 1:length(sessions)
%         session = sessions(sCounter);
%         for counter = 1:session.SessionData.nTrials
%             TE.filename{tcounter} = session.filename;
%             TE.trialNumber(tcounter) = counter;
%             TE.trialType(tcounter) = session.SessionData.TrialTypes(counter);
%             TE.trialOutcome(tcounter) = session.SessionData.TrialOutcome(counter);
%             TE.epoch(tcounter) = session.SessionData.Epoch(counter);
%             tcounter = tcounter + 1; % don't forget :)            
%         end
%     end