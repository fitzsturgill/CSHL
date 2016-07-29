function timerTrigger
	global state
	state.cycle.cycleStatus=2; 	% acquiring
	state.cycle.lastPositionUsed = state.cycle.currentCyclePosition;
	
	try
		disp(['Triggered at ' clockToString(clock) ', ' num2str(etime(clock,state.internal.triggerTime)) ' seconds since last acquisition.']);
	catch
		disp(['Triggered at ' clockToString(clock)]);
	end

	dioTrigger;
    
    if state.cycle.mcOlfOn % if olfactometer is actuated during this cycle
        mcOlfTrigger;
    end