function mcAcqUpdateChannelNames            
    global state
    
        channelOrder = state.phys.mcAcq.mcChannelOrder;
            
        for i = 1:state.phys.mcAcq.totalChannels
            if i<=state.phys.mcAcq.mcNChannels
                state.phys.mcAcq.channel(i).Name = ['mc' num2str(channelOrder(i))];
                state.phys.mcAcq.mcChannelNames{1, i} = ['mc' num2str(channelOrder(i))];
                state.phys.mcAcq.channel(i).Show = 1;
            else
                state.phys.mcAcq.channel(i).Name = ['aux' num2str(i)];
                state.phys.mcAcq.mcChannelNames{1, i} = ['aux' num2str(i)];
            end            
        end
