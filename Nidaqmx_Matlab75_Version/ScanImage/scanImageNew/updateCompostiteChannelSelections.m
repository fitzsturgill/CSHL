function updateCompostiteChannelSelections
	global state
	
	state.internal.compositeChannelSelections=zeros(1,10+state.init.maximumNumberOfInputChannels);
	
	selectedChannels=[state.internal.redCompositeChannel state.internal.greenCompositeChannel state.internal.blueCompositeChannel];
	
	for channelCounter=1:3
		if selectedChannels(channelCounter)>1
			if selectedChannels(channelCounter)<=state.init.maximumNumberOfInputChannels+1
				state.internal.compositeChannelSelections(selectedChannels(channelCounter)-1)=channelCounter;
			elseif selectedChannels(channelCounter)>state.init.maximumNumberOfInputChannels
				state.internal.compositeChannelSelections(10+selectedChannels(channelCounter)-state.init.maximumNumberOfInputChannels-1)=channelCounter;
			end
		else
			state.acq.compositeData(:,:,channelCounter)=0;
		end
	end