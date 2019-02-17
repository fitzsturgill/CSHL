function [avg avgX avgSEM allTrials ] = sessionAvg(SessionData, type, outcome, linespec, ax)
% make SessionData global eventually???

    
    sample_rate = 10000;
    fCh = 1;
    baselineDuration = 2; % duration of baseline period prior to start of C.S.
    trials = fsFilterTrials(SessionData, type, outcome);
    nTrials = length(trials);
    periStimF = cell(nTrials, 1);
    startIndices=[];
    for i = 1:length(trials)
        trial = trials(i);
        %             % GCaMP response
            stimulus_start_time = SessionData.RawEvents.Trial{trial}.States.DeliverStimulus(1);
            startIndex = round((stimulus_start_time - 2) * sample_rate); %start 2 seconds prior to stim start
%             sample_indices = round(stimulus_start_time * sample_rate + psth_flanking_samples(:) );
%             available_samples = size( SessionData.NidaqData{trial} , 1) - stimulus_start_time * sample_rate;
%             if available_samples < 0,
%                 fprintf('available_samples < 0 for some reason (%.0f) :: trial %.0f :: stimulus\n',available_samples,trial)
%                 continue
%             end
%             end_sample = round(min(length(psth_flanking_samples),available_samples));
            % Symmetric window around start
            %stimulus_zeroing_window = (stimulus_start_time * sample_rate) + round(linspace(-floor(zeroing_samples/2),floor(zeroing_samples/2)-1,zeroing_samples));
            % Causal window around start
%             stimulus_zeroing_window = (stimulus_start_time * sample_rate) + round(linspace(-floor(sample_rate),-1,sample_rate));
            % Populate peri_stim_GCaMP matrix with response normalized by mean around t=0;

            trialData = SessionData.NidaqData{trial}(startIndex:end , fCh )'; % make it a row vector
%             trialData = SessionData.NidaqData{trial}(: , fCh )'; % make it a row vector
            if isempty(trialData)
                disp(['No Data in trial #' num2str(trial)]);
            end
            periStimF{i} = trialData;
    end
    
    [avg avgSEM allTrials] = avgTrialsF(periStimF);
    avgX = (1:length(avg)) .* (1/sample_rate);
    avgX = avgX - 2;
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
        allTrials = allTrials - baselines;  % baseline the trials
        avg = mean(allTrials);
        avgN = nTrials;
        avgSEM = std(allTrials) ./ sqrt(avgN);
    end

end
    

 