function mcAcqUpdateChannelPlots
    global state gh  mcData
    
%     if state.mcViewer.showThresh
%         show = 'on';
%     else
%         show = 'off';
%     end
        
    validChannels=find(mcChannelFieldToVector('Show'));
    state.phys.mcAcq.displayedChannels = validChannels;
    nValidChannels = length(validChannels);
    
    
% filter multiChannel Data...
% I'm filtering each channel individually to allow for different filter
% settings for different channels (even though this is likely slower)
    
%      lowPass = mcChannelFieldToVector('LowPass');
%      highPass = mcChannelFieldToVector('HighPass');
%      showFilter = mcChannelFieldToVector('ShowFilter');
     
    for i = 1:nValidChannels
        channel = validChannels(i);
        if state.phys.mcAcq.channel(channel).ShowFilter
            % no need to refilter data if channel filtering parameters match
            % global filtering parameters
            if channel <= state.phys.mcAcq.mcNChannels && ((state.phys.mcAcq.channel(channel).LowPass == 0 && state.phys.mcAcq.channel(channel).HighPass == 0) || (state.phys.mcAcq.channel(channel).LowPass == state.phys.mcAcq.globalLowPass && state.phys.mcAcq.channel(channel).HighPass == state.phys.mcAcq.globalHighPass))
                state.phys.mcAcq.displayData(:, channel) = state.phys.mcAcq.filteredData(:, channel);
            else
                state.phys.mcAcq.displayData(:, channel) = mcFilter(mcData.data(:, channel),...
                state.phys.mcAcq.channel(channel).LowPass, state.phys.mcAcq.channel(channel).HighPass, state.phys.mcAcq.mcInputRate);
            end
        end
    end
    
    
% Plot multichannel channels selected for display
    for i = 1:nValidChannels
        channel = validChannels(i);
        set(state.phys.mcAcq.lines(1, i),...
            'Ydata', state.phys.mcAcq.displayData(:, channel),...
            'Xdata', state.phys.mcAcq.displayXData...
            );
        
    end