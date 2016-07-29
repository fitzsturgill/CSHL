function closeChannelGUI
	global state

	if state.internal.channelChanged
		hideGUI('gh.channelGUI.figure1');
		updateChannelFlags;
		applyImagingInputParams;	
	else
		hideGUI('gh.channelGUI.figure1');
		state.internal.channelChanged=0;
	end
	