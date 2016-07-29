function setPcellsToDefault
	global state
	
	if ~state.pcell.pcellOn
		return
	end
	
	stop([state.daq.pcellFocusOutput state.daq.pcellGrabOutput]);
	vec=[];

	for counter=1:state.pcell.numberOfPcells
		vec(counter)=powerToPcellVoltage(getfield(state.pcell, ['pcellDefaultLevel' num2str(counter)]), counter);
		vec(counter+state.pcell.numberOfPcells) = 5 * state.shutter.closed;
	end

	putsample(state.daq.pcellFocusOutput, vec);	


			