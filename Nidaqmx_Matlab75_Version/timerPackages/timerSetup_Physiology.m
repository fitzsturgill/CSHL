function timerSetup_Physiology
	global state

	if state.phys.internal.abort
%		abortPhysiology;
		return
	end
	
	try
		readTelegraphs;
		updateMinInCell;
	catch
		disp(['timerSetup_Physiology: ' lasterr]);
	end
			
	setupPhysDaqPulse;
    mcOlfSetup;
