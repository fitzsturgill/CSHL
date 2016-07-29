function flushPhysData
	global state
	
	try
		if get(state.phys.daq.outputDevice, 'SamplesAvailable')>0
			start(state.phys.daq.outputDevice);
			stop(state.phys.daq.outputDevice);
	
			while ~strcmp(get(state.phys.daq.outputDevice, 'Running'), 'Off')
				pause(.001);
			end
		end
	catch
		disp(['flushPhysData: ' lasterr]);
	end
	
% 	try
% 			
% 		start(state.phys.daq.outputDevice);
% 		stop(state.phys.daq.outputDevice);
% 	
% 		while ~strcmp(state.phys.daq.outputDevice.Running, 'Off');
% 			pause(.001);
% 		end
% 			catch
% 		disp('flushPhysData: Error.  Assume aborted.');
% 	end
