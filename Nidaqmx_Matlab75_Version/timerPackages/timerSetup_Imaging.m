function timerSetup_Imaging
	global state
	applyChangesToOutput;			
	
	if state.init.autoReadPMTOffsets
		getPMTOffsets;
	end
	mp285Flush;
	if state.acq.numberOfZSlices > 1	
		state.internal.initialMotorPosition=updateMotorPosition;
	else
		state.internal.initialMotorPosition=[];
	end		
	turnOffMenus;
	
	
		
