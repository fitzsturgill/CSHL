function flushPhysAuxData
	global state
	disp ' called flushphysauxdaat'
	try
			
		start(state.phys.daq.auxOutputDevice);
		stop(state.phys.daq.auxOutputDevice);
	
		while ~strcmp(state.phys.daq.auxOutputDevice.Running, 'Off');
			pause(.001);
		end
	catch
		disp('flushPhysAuxData: Error.  Assume aborted.');
	end
