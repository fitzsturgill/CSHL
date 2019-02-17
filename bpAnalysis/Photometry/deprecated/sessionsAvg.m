function [avg avgX avgSEM allTrials ] = sessionsAvg(sessions, type, outcome, linespec, ax)
% make sessions global eventually???
% kludgy function to average over sessions for figure for research proposal
    
    sample_rate = 10000;
    fCh = 1;
    baselineDuration = 2; % duration of baseline period prior to start of C.S.
    
    periStimF = {};
    for j=1:length(sessions)
        session = sessions(j).SessionData;
        trials = fsFilterTrials(session, type, outcome);


        for i = 1:length(trials)
            trial = trials(i);
            %             % GCaMP response
                stimulus_start_time = session.RawEvents.Trial{trial}.States.DeliverStimulus(1);
                startIndex = round((stimulus_start_time - 2) * sample_rate); %start 2 seconds prior to stim start


                trialData = session.NidaqData{trial}(startIndex:end , fCh )'; % make it a row vector
                trialData = smooth(trialData', 5000)';
                
                if isempty(trialData) || numel(trialData) < (7 * sample_rate)
                    continue
                else
                    periStimF{end + 1, 1} = trialData;
                end
        end
    end

        
        
        
        
        [avg avgSEM allTrials] = avgTrialsF(periStimF);
        avgX = (1:length(avg)) .* (1/sample_rate);
        avgX = avgX - 2;
        
        % I need to downsample
        avgX = decimate(avgX, 100);
        avg = decimate(avg, 100);
        avgSEM = decimate(avgSEM, 100);
    %    plot(xData, allTrials'); hold on

        boundedline(avgX, avg, avgSEM, linespec, ax, 'alpha');
    


    function [avg avgSEM allTrials] = avgTrialsF(trialsF)
        % baseline samples currently coded as 2sec * 10000Hz = 20000
        baselineSamples = baselineDuration * sample_rate;
        nTrials = size(trialsF, 1);
        minSamples = min(cellfun(@length, trialsF));
        %initialize
        allTrials = zeros(nTrials, minSamples);

        for i=1:nTrials
            allTrials(i, :) = trialsF{i}(1:minSamples); 
        end

        baselines = mean(allTrials(:, 1:baselineSamples), 2);
        baselines = repmat(baselines, 1, size(allTrials, 2));
        allTrials = allTrials - baselines;  % baseline the trials, deltaF
        allTrials = allTrials ./ baselines; % normalize, the trials, deltaF/F
        avg = mean(allTrials);
        avgN = nTrials;
        avgSEM = std(allTrials) ./ sqrt(avgN);
    end

end
    

 