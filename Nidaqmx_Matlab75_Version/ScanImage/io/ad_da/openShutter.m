function openShutter
	global state

	putvalue(state.daq.shutterLine, state.shutter.open);