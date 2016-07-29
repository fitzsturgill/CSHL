function out=timerFirstSetup_Physiology
	out=0;
	global state gh
	
	state.phys.internal.abort=0;
	state.phys.internal.first=1;
	state.phys.scope.changedScope=0;
	state.phys.internal.outputNeedsUpdate=1;
	
	
