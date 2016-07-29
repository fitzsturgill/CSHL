function putDataFocus
	global state
	
	putdata(state.daq.focusOutput, state.acq.rotatedMirrorData);	
	putdata(state.daq.pcellFocusOutput, state.acq.pcellPowerOutput);		
