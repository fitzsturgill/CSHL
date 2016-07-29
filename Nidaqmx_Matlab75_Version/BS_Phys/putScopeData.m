function putScopeData
	global state
	unique(state.phys.scope.output);
	putdata(state.phys.daq.scopeOutputDevice, state.phys.scope.output)
	