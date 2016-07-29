function abortFocus
	global state gh
	state.internal.abortActionFunctions=0;
	
	stop([state.daq.pcellFocusOutput state.daq.focusOutput state.daq.focusInput]);

	while ~strcmp(get(state.daq.pcellFocusOutput, 'Running'), 'Off')
		pause(0.001);
	end
	
	setPcellsToDefault;
	setStatusString('Aborting Focus...');
	
	closeShutter;
	set(gh.standardModeGUI.focusButton, 'Enable', 'off');

	while ~any(strcmp(get([state.daq.pcellFocusOutput state.daq.focusOutput state.daq.focusInput], 'Running'), 'Off'))
		pause(0.001);
	end
	
	setPcellsToDefault;

	flushFocusData;

	if get(state.daq.focusInput, 'SamplesAvailable')>0
		try
			flushdata(state.daq.focusInput);
		catch
			disp('abortFocus: error in input flush data.  proceeding...');
		end
	end
	
	mp285Flush;

	
	set(gh.standardModeGUI.focusButton, 'String', 'FOCUS');
	set(gh.standardModeGUI.focusButton, 'Enable', 'on');
%	set(gh.standardModeGUI.startLoopButton, 'Visible', 'On');

	if ~state.internal.looping
		set(gh.standardModeGUI.grabOneButton, 'Visible', 'On');
		turnOnMenus;
		state.internal.status=0;
		applyChangesToOutput;
	else
		state.internal.status=4;
		applyChangesToOutput;
		turnOffMenus;

		resetCounters;
		state.internal.abortActionFunctions=0;
		setStatusString('Resuming loop...');
	
		updateGUIByGlobal('state.internal.frameCounter');
		updateGUIByGlobal('state.internal.zSliceCounter');

		state.internal.abort=0;
		state.internal.currentMode=3;

		mainLoop;
	end

	setStatusString('');
	timerReleasePause('Imaging');
	
