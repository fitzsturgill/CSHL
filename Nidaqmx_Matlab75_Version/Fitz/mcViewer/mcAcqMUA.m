function mcAcqMUA
    global state mcData 
   


    
    % for each epoch an average is created if one or more channels have MUA
    % selected.  Copy all functionality is disabled so that you don't
    % inadvertantly reset averages
    
    % the averages are in wave format, i.e. part of the WAVE class created
    % by BLS emulating the wave functionality in IGOR
    % A grand average MUA is also created for which the trial MUAs of the
    % selected channels are first averaged and then averaged into the
    % grand average
    
    % averages are grouped according to cycle position
    
    % thresholds are set individually based upon the SD of the first trace
    % acquired after beginning a new epoch or resetting the average
    
    % if you alter the channels to be included in the MUA then the MUA is
    % reset
    
    % see additional functions for displaying the waves.  Wave plots are
    % automatically updated
    

        
    includeChannels = find(mcChannelFieldToVector('MUAInclude'));
    if isempty(includeChannels)
        return
    end
   binSize = 500; % 200 msec
   binSize = binSize/1000; % convert to seconds
    binEdges = 0:binSize:max(state.phys.mcAcq.displayXData);
%     iterate over channels

    % determine if all channels average exists or is empty
    wnAll = MUAAvgName;  %'All'
    if ~iswave(wnAll) || state.phys.mcAcq.MUA.resetMUA
        newAverage(wnAll)
        setwave(wnAll, 'xscale', [0 binSize]);
    end
    

    for i = 1:length(includeChannels)
    
    % first determine if average already exists or is empty, then set
    % automatic thresholds only if average does not prexist or is empty
        channel = includeChannels(i);
        wn = MUAAvgName(channel);
        if ~iswave(wn) || state.phys.mcAcq.MUA.resetMUA
            if ~iswave(wn)
                newAverage(wn);
            end
            setwave(wn, 'xscale', [0 binSize]);
            % determine threshold, negative going
            threshold = state.phys.mcAcq.MUA.means(1, channel) - state.phys.mcAcq.MUA.SD(1, channel) * 3;
            state.phys.mcAcq.MUA.thresholds(1, channel) = threshold;   % later I may implement manual threshold mode
            state.phys.mcAcq.MUA.autoThresholds(1, channel) = threshold;
        end
        
        % determine spike times
        data = state.phys.mcAcq.filteredData(:, channel)'; % convert to row vector
        crossDownIndices = diff((data - state.phys.mcAcq.MUA.means(1, channel)) < state.phys.mcAcq.MUA.thresholds(1, channel)) == -1 ; 
        crossDown = state.phys.mcAcq.displayXData(crossDownIndices)';
        state.phys.mcAcq.MUA.spikeTimes = crossDown;
        
        % make trial histogram
        spikeHist = histc(state.phys.mcAcq.MUA.spikeTimes, binEdges);
        spikeHist = spikeHist / binSize;% convert to Hz
        % average together channels for All histogram
        if i == 1
            avgAllData = zeros(size(spikeHist));
        end
        avgAllData = avgAllData + spikeHist;
        % avgin
        waveo('tempHistWave', spikeHist) % convert to wave
        setwave('tempHistWave', 'xscale', [0 binSize]);
        avgin('tempHistWave', wn);
    end
    
    if state.phys.mcAcq.MUA.resetMUA
        state.phys.mcAcq.MUA.resetMUA=0;
        updateGUIByGlobal('state.phys.mcAcq.MUA.resetMUA'); % turn off resetMUA checkbox
    end
    
    avgAllData = avgAllData ./ length(includeChannels);
    waveo('tempHistAllWave', avgAllData);
    setwave('tempHistAllWave', 'xscale', [0 binSize]);
    avgin('tempHistAllWave', wnAll);
        
        
        
    % make a wave to contain histogram for current acquisition-  i.e.
    % MC_hist
    % then use avgin to create the avg, e.g.  avgin(MC_hist,
    % 'PirC_12_c1_avg');
    % components stored in average will be useless as I'm not going to fill
    % up memory with waves corresponding to each acquisition
    % however, avgin, resetavg and plot-related functionality will be
    % useful....
    
    
    
    
    

