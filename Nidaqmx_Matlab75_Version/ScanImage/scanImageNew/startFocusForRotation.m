function startFocusForRotation
	global state
	if get(state.daq.focusInput, 'SamplesAvailable')>0 
		try
			flushData(state.daq.focusInput);
		catch
			disp(['startFocusForRotation: ' lasterr]);
		end
	end

	putDataFocus;
	
	resetCounters;
		
	startFocus;
	openShutter;
	state.internal.abortActionFunctions=0;
	
	diotrigger;