function executeStartLoopCallback;

	global state gh
	state.internal.whatToDo=3;


	val=get(gh.mainControls.startLoopButton, 'String');
	state.internal.cyclePaused=0;
	
	if strcmp(val, 'LOOP')
		if strcmp(get(gh.basicConfigurationGUI.figure1, 'Visible'), 'on')
			beep;
			setStatusString('Close ConfigurationGUI');
			return
		end
			
		if ~savingInfoIsOK
			return
		end
		
		mp285Flush;
		set(gh.mainControls.startLoopButton, 'String', 'ABORT');
		set(gh.mainControls.grabOneButton, 'Visible', 'Off');
		turnOffMenus;
		
		resetCounters;
		state.internal.abortActionFunctions=0;

		setStatusString('Starting loop...');

		state.internal.firstTimeThroughLoop=1;
		state.acqParams.triggerTime=clock;
		state.internal.abort=0;
		if state.timer.timerActive
			setStatusString('Setting up packages...');
			timerCallPackageFunctions('FirstSetup');
		end
		mainLoop;
	else
		set(gh.mainControls.startLoopButton, 'Enable', 'off');
		state.internal.looping=0;
		state.internal.abort=1;
		abortGrab;
		if state.timer.timerActive
			timerCallPackageFunctions('Abort');
		end
		
		setStatusString('Stopping loop...');

		set([gh.mainControls.focusButton gh.mainControls.grabOneButton], 'Visible', 'On');
		turnOnMenus;
		set(gh.mainControls.startLoopButton, 'String', 'LOOP');
		set(gh.mainControls.startLoopButton, 'Enable', 'on');
		setStatusString('');
	end



