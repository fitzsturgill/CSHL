function setTimerStatusString(st)

	global state
	state.timer.statusString=st;
	updateGUIByGLobal('state.timer.statusString');
	
	