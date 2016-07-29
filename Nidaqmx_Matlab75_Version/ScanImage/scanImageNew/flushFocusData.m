function flushFocusData
	global state
	
	try
		deviceList=[];
		if state.pcell.pcellOn
			if get(state.daq.pcellFocusOutput, 'SamplesAvailable')>0
				deviceList=state.daq.pcellFocusOutput;
			end
		end
		if get(state.daq.focusOutput, 'SamplesAvailable')>0
			if isempty(deviceList)
				deviceList=state.daq.focusOutput;
			else
				deviceList(end+1)=state.daq.focusOutput;
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
		disp(['flushFocusData: ' lasterr]);
	end

