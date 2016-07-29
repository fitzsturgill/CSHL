function abortPCells
	global state
	if state.pcell.pcellOn
		stop([state.daq.pcellFocusOutput state.daq.pcellGrabOutput]);
	end
	setPcellsToDefault;
