function states = bpAddStatesAsTrialEvents(session)
    
    stateList = {};
    nTrials = session.SessionData.nTrials;
    rawEvents = session.SessionData.RawEvents.Trial; % cell array
    states = struct();
    for trial = 1:nTrials
        trialFields = fieldnames(rawEvents{trial}.States);
        for i = 1:length(trialFields)
            field = trialFields{i};
            if ~ismember(field, stateList)
                states.(field) = repmat({NaN}, nTrials, 1);
                stateList{end + 1} = field;
            end
            states.(field){trial} = rawEvents{trial}.States.(field);
        end
    end
            
        
        