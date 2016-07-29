function timerAbort_Imaging
	global state gh

	state.internal.abortActionFunctions=1;
	state.internal.status=0;
	
	if state.acq.externalTrigger
		if strcmp(get(state.daq.grabInput, 'Running'),'Off') 
			abortGrab;
		elseif get(state.daq.grabInput, 'TriggersExecuted')>0
			state.internal.abortActionFunctions=1;
		else
			abortGrab;
		end
	else
		if strcmp(get(state.daq.grabInput, 'Running'),'Off') | get(state.daq.grabInput, 'TriggersExecuted')==0
			abortGrab;
		else
			state.internal.abortActionFunctions=1;
		end
	end
	
