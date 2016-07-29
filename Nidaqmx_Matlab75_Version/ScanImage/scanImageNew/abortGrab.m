function abortgrab
	global state gh

	stop([state.daq.pcellGrabOutput state.daq.grabOutput state.daq.grabInput]);

	while ~strcmp(get(state.daq.pcellGrabOutput, 'Running'), 'Off')
		pause(0.001);
	end
	state.internal.abortActionFunctions=0;
	setPcellsToDefault;

	setStatusString('Aborting Grab...');

	set(gh.standardModeGUI.grabOneButton, 'Enable', 'off');

	while ~any(strcmp(get([state.daq.pcellGrabOutput state.daq.grabOutput state.daq.grabInput], 'Running'), 'Off'))
		pause(0.001);
	end

	state.internal.status=0;
	try
		flushData(state.daq.grabInput);
	catch
%		disp('abortGrab: error in input flush data.  proceeding...');
	end
		
	mp285Flush;
	
	executeGoHome;

	set(gh.standardModeGUI.grabOneButton, 'String', 'GRAB');
	set(gh.standardModeGUI.grabOneButton, 'Enable', 'on');
	set(gh.standardModeGUI.focusButton, 'Visible', 'On');

	timerSetPackageStatus(0, 'Imaging');
	timerCheckIfAllAborted;
	resetCounters;
	setStatusString('');
	
