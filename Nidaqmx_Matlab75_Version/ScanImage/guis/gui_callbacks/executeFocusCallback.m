function executeFocusCallback
	global state gh
	
	if strcmp(get(gh.standardModeGUI.focusButton, 'String'), 'FOCUS')
		if strcmp(get(gh.basicConfigurationGUI.figure1, 'Visible'), 'on')
			beep;
			setStatusString('Close ConfigurationGUI');
			return
		end
		
		if timerGetPackageStatus('Imaging')
			beep;
			setStatusString('Running. Stop processes');
			return
		end

		timerRequestPause('Imaging');
		setStatusString('Focusing...');
		set(gh.standardModeGUI.focusButton, 'String', 'ABORT');

		set(gh.standardModeGUI.grabOneButton, 'Visible', 'Off');

		if state.internal.looping
			state.internal.cyclePaused=1;
		end
		turnOffMenus;
		
		if state.init.autoReadPMTOffsets
			getPMTOffsets;
		end
		
 		if state.internal.updatedZoomOrRot % need to reput the data with the approprite rotation and zoom.
			if get(state.daq.focusOutput, 'SamplesAvailable')>0
				flushFocusData;
			end
			if get(state.daq.grabOutput, 'SamplesAvailable')>0
				flushGrabData;
			end
		
			putDataFocus;
			putDataGrab;
			state.internal.updatedZoomOrRot=0;
 		end
		
		mp285Flush;
		resetCounters;
		
		state.internal.abortActionFunctions=0;

		startFocus;
		openShutter;
		diotrigger(0);
	else
		if strcmp(get(state.daq.focusInput, 'Running'),'Off')
			abortFocus;
		else
			state.internal.abortActionFunctions=1;
		end
	end
	

