function setPhysStatusString(status)
	global state
	state.phys.internal.statusString=status;
	updateGuiByGlobal('state.phys.internal.statusString');
	