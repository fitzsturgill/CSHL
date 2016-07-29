function flushGrabData
	global state
	
	try
		deviceList=[];
		if state.pcell.pcellOn
			if get(state.daq.pcellGrabOutput, 'SamplesAvailable')>0
				deviceList=state.daq.pcellGrabOutput;
			end
		end
		if get(state.daq.grabOutput, 'SamplesAvailable')>0
			if isempty(deviceList)
				deviceList=state.daq.grabOutput;
			else
				deviceList(end+1)=state.daq.grabOutput;
			end
		end
		
		if ~isempty(deviceList)
			start(deviceList);
			stop(deviceList);
	
			while ~any(strcmp(get(deviceList, 'Running'), repmat('Off', length(deviceList), 1)))
				pause(.001);
			end
		end
	catch
		disp(['flushGrabData: ' lasterr]);
	end
