function timerStart_Imaging
	global state gh

	timerSetPackageStatus(1, 'Imaging');
	resetCounters;
	state.internal.abortActionFunctions=0;
		
	set(gh.standardModeGUI.grabOneButton, 'String', 'ABORT');
	set(gh.standardModeGUI.focusButton, 'Visible', 'Off');

	state.files.lastAcquisition=state.files.fileCounter;

	deviceList=[state.daq.grabInput state.daq.grabOutput];
	if state.pcell.pcellOn
		deviceList(end+1)=state.daq.pcellGrabOutput;
	end
	start(deviceList)
	state.internal.status=3;
	state.internal.lastTaskDone=3;




 	
