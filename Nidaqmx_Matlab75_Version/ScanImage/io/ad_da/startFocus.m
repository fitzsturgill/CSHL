function startFocus
	global state
	
	state.internal.status=2;
	state.internal.lastTaskDone=2;
	if state.pcell.pcellOn
		start([state.daq.focusInput state.daq.focusOutput state.daq.pcellFocusOutput]);
	else
		start([state.daq.focusInput state.daq.focusOutput]);
	end



 