function abort=timerWaitWhilePaused
	abort=0;
	global state
	while timerPausedStatus
		timerPausedStatus
		setStatusString('LOOP PAUSED');
		pause(0.1);
		if state.timer.abort
			abort=1;
			return
		end	
	end	
	if state.timer.abort
		abort=1;
	end	
