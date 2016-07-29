function startGrab
	global state

	deviceList=[state.daq.grabInput state.daq.grabOutput];
	if state.pcell.pcellOn
		deviceList(end+1)=state.daq.pcellGrabOutput;
	end
	start(deviceList)
	state.internal.status=3;
	state.internal.lastTaskDone=3;
