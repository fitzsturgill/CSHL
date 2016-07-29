function stopAndRestartFocus
	global state

	state.internal.pauseAndRotate=0;

	deviceList=[state.daq.focusInput state.daq.focusOutput];
	if state.pcell.pcellOn
		deviceList(end+1)=state.daq.pcellFocusOutput;
	end
	
	stop(deviceList);
	while ~any(strcmp(get(deviceList, 'Running'), 'Off'))
		pause(0.001);
	end	

	for counter=1:state.pcell.numberOfPcells
		vec(counter)=powerToPcellVoltage(getfield(state.pcell, ['pcellDefaultLevel' num2str(counter)]), counter);
		vec(counter+state.pcell.numberOfPcells) = 5 * state.shutter.closed;
	end
	putsample(state.daq.pcellFocusOutput, vec);	

	flushFocusData;

	if get(state.daq.focusInput, 'SamplesAvailable')>0
		try
			flushData(state.daq.focusInput);
		catch
			disp(lasterr);
		end
	end

	putDataFocus;
	
	state.internal.stripeCounter=0;
	state.internal.focusFrameCounter = 1;
			
	startFocus;
	diotrigger;

	