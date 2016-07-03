function [sessions, TE, allTE] = makeTE_SO_RewardPunish_odor(sessions, prepRaster)
    
    if nargin < 2
        prepRaster = 0; % kludge, if set to 1, add phRaster field to 
%         allTE, this should only be set if trials are of constant length
%         (for now)
    end
    
    TE = struct();   
    %%
    for pos = 1:length(sessions);
%         sessions = bpLoadSession(sessionNames{pos}, sessionPath);
        sessions(pos).SessionData = demodulateSession(sessions(pos).SessionData);
        nTrials = sessions(pos).SessionData.nTrials;
        TE(pos).Trial = 1:length(sessions(pos).SessionData.TrialTypes);
        TE(pos).TrialTypes = sessions(pos).SessionData.TrialTypes;
        TE(pos).TrialOutcome = sessions(pos).SessionData.TrialOutcome;
        TE(pos).blLicks = bpCountEventByStates(sessions(pos), 'Port1In', 'Cue', 'window', [-2 0]); % 2 seconds prior to cue
        TE(pos).anticipatoryLicks1 =  bpCountEventByStates(sessions(pos), 'Port1In', 'Cue', 'endState', 'Delay'); % all licks cue + delay
        TE(pos).anticipatoryLicks2 = bpCountEventByStates(sessions(pos), 'Port1In','Delay', 'referenceFromEnd', 1, 'window', [-2 0]);    % fixed 2 second window prceding US
        TE(pos).cueLicks = bpCountEventByStates(sessions(pos), 'Port1In', 'Cue', 'window', [0.4 1.4]); % cue licks, matches phCuePeak
        TE(pos).cueLicksLong = bpCountEventByStates(sessions(pos), 'Port1In', 'Cue', 'window', [0.4 2.4]); % cue licks, matches phCuePeak    
        TE(pos).delayLicks = bpCountEventByStates(sessions(pos), 'Port1In', 'Delay', 'window', [0 1]); % delay licks (just 1st second), matches phDelayLicks, note that delay is somewhat variable
        TE(pos).usLicks = bpCountEventByStates(sessions(pos), 'Port1In', 'Delay', 'referenceFromEnd', 1, 'window', [0 1]); % us licks, matches phUsPeak
        TE(pos).Photometry = processTrialAnalysis_Photometry(sessions(pos));
        TE(pos).states = bpAddStatesAsTrialEvents(sessions(pos));
        if isfield(sessions(pos).SessionData, 'Epoch');
            TE(pos).epoch = sessions(pos).SessionData.Epoch;
        else
            TE(pos).epoch = ones(1, nTrials);
        end
        TE(pos).fileName = repmat({sessions(pos).filename}, nTrials, 1);
        TE(pos).phCuePeak = bpCalcPeak_dFF(TE(pos).Photometry, 1, [0.4 1.4], TE(pos).states.Cue); % cue is always 1s
        TE(pos).phCuePeakLong = bpCalcPeak_dFF(TE(pos).Photometry, 1, [0.4 2.4], TE(pos).states.Cue); % cue is always 1s
        TE(pos).phDelayPeak = bpCalcPeak_dFF(TE(pos).Photometry, 1, [0 1], TE(pos).states.Delay); % my delay is variable (1s initially, 2s for most), still I'm just taking 1st second
        TE(pos).phUsPeak = bpCalcPeak_dFF(TE(pos).Photometry, 1, [0 1], TE(pos).states.Delay, 'referenceFromEnd', 1);
        TE(pos).phBaseline = bpCalcPeak_dFF(TE(pos).Photometry, 1, [-2 0], TE(pos).states.Cue);    % baseline 2s before cue
    end
    
    
    allTE = struct(...
        'Trial', [],...
        'TrialTypes', [],...
        'TrialOutcome', [],...
        'blLicks', [],...
        'anticipatoryLicks1', [],...
        'anticipatoryLicks2', [],...
        'cueLicks', [],...
        'cueLicksLong', [],...
        'delayLicks', [],...
        'usLicks', [],...
        'epoch', [],...
        'fileName', {},...
        'phCuePeak', [],...
        'phCuePeakLong', [],...
        'phDelayPeak', [],...
        'phUsPeak', [],...
        'phBaseline', []...
        );
        
    if prepRaster % kludgy
        sessionsCount = zeros(1, length(TE));
        samplesCount = zeros(1, length(TE));
        for counter = 1:length(TE)
            sessionsCount(counter) = TE(counter).Trial(end);
            samplesCount(counter) = size(TE(counter).Photometry.data(1).dFF, 2);
        end
        sessionsCount = sum(sessionsCount);
        samplesCount = max(samplesCount);
        allTE(1).phRaster = zeros(sessionsCount, samplesCount); % initialize
    end
    %     'dFF', {},... % ch1
    nextSessionIndex = 1;
    % WTF am I doing here with the allTE(1).TrialTypes = [allTE.TrialTypes;
    for counter = 1:length(TE)
        nTrials = TE(counter).Trial(end);
        allTE(1).Trial = [allTE(1).Trial; (1:nTrials)'];
        allTE(1).TrialTypes = [allTE.TrialTypes ; TE(counter).TrialTypes(1:length(TE(counter).TrialOutcome))'];
        allTE(1).TrialOutcome = [allTE.TrialOutcome; TE(counter).TrialOutcome'];
        allTE(1).blLicks = [allTE.blLicks; TE(counter).blLicks.rate'];
        allTE(1).anticipatoryLicks1 = [allTE.anticipatoryLicks1; TE(counter).anticipatoryLicks1.rate'];
        allTE(1).anticipatoryLicks2 = [allTE.anticipatoryLicks2; TE(counter).anticipatoryLicks2.rate'];
        allTE(1).cueLicks = [allTE.cueLicks; TE(counter).cueLicks.rate'];
        allTE(1).cueLicksLong = [allTE.cueLicksLong; TE(counter).cueLicksLong.rate'];    
        allTE(1).delayLicks = [allTE.delayLicks; TE(counter).delayLicks.rate'];
        allTE(1).usLicks = [allTE.usLicks; TE(counter).usLicks.rate'];
        allTE(1).epoch = [allTE.epoch ; TE(counter).epoch'];
        allTE(1).fileName = [allTE.fileName ; TE(counter).fileName];
        allTE(1).phCuePeak = [allTE.phCuePeak ; TE(counter).phCuePeak.data];
        allTE(1).phCuePeakLong = [allTE.phCuePeakLong ; TE(counter).phCuePeakLong.data];    
        allTE(1).phDelayPeak = [allTE.phDelayPeak ; TE(counter).phDelayPeak.data];    
        allTE(1).phUsPeak = [allTE.phUsPeak ; TE(counter).phUsPeak.data];  
        allTE(1).phBaseline = [allTE.phBaseline ; TE(counter).phBaseline.data];
        if prepRaster
            endPoint = size(TE(counter).Photometry.data(1).dFF, 2);
            allTE(1).phRaster(nextSessionIndex:nextSessionIndex + nTrials - 1, 1:endPoint) = TE(counter).Photometry.data(1).dFF;
            nextSessionIndex = nextSessionIndex + nTrials;
        end
        
    end
    %%
    % convert fileName into indices
    allNames = unique(allTE.fileName);
    sessionIndex = zeros(size(allTE.TrialOutcome)); % now unique names are indices to seprates sessions
    for counter = 1:length(allNames)
        sessionIndex(strcmp(allNames{counter}, allTE.fileName)) = counter;
    end
    allTE.sessionIndex = sessionIndex;
    allTE.trialIndex = 1:length(allTE.epoch);
