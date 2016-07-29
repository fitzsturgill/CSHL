function timerWait_Physiology
	global state
	
	if timerGetPackageStatus('Physiology')
		return
	end
	
	try
		state.phys.internal.timer=state.internal.secondsCounter;
		updateGuiByGlobal('state.phys.internal.timer');
		readTelegraphs;
		if state.phys.scope.changedScope
			state.phys.scope.changedScope=0;
			setupCyclePosition;
			setUpPhysDaqPulse;						
			readBaseline;
		else
			readBaseline;
		end
		updateMinInCell;
    catch
		disp(['timerWait_Physiology: ' lasterr]);
    end