function abortPhysiology
	global state gh

	state.phys.internal.abort=1;

	inputRunning=0;
	if ~strcmp(state.phys.daq.inputDevice.Running, 'Off')
		stop(state.phys.daq.inputDevice);
		inputRunning=1;
	end

	outputRunning=0;
	if ~strcmp(state.phys.daq.outputDevice.Running, 'Off')
		outputRunning=1;
		stop(state.phys.daq.outputDevice);
	end

	if state.phys.daq.auxOutputBoardIndex & any(state.cycle.lastUsedAuxPulses)
		if ~strcmp(state.phys.daq.auxOutputDevice.Running, 'Off')
			auxOutputRunning=1;
			stop(state.phys.daq.auxOutputDevice);
		end
		while ~strcmp(state.phys.daq.auxOutputDevice.Running, 'Off')
			pause(0.001);
		end	
		if get(state.phys.daq.auxOutputDevice, 'SamplesAvailable')>0
			flushPhysAuxData;
		end
	end
		
	while ~strcmp([state.phys.daq.inputDevice.Running state.phys.daq.outputDevice.Running], ...
		['Off' 'Off']);
		pause(0.001);
	end	

	if get(state.phys.daq.outputDevice, 'SamplesAvailable')>0
		flushPhysData;
	end

	if inputRunning
		flushData(state.phys.daq.inputDevice);
	end

	set(gh.physControls.startButton, 'String', 'START');
	pause(0.001);
	set(gh.physControls.startButton, 'Enable', 'on');
	set(gh.scope.start, 'Enable', 'on');
    
	setPhysStatusString('Ready');	
	timerSetPackageStatus(0, 'Physiology');
	timerCheckIfAllAborted;
    
    % FS MOD
    
    % make sure that mcAcquisition is switched off
    if state.phys.mcAcq.mcInputBoardIndex
        stop(state.phys.mcAcq.mcInputDevice);
    end
    % make sure that odor is switched off
    if state.phys.mcAcq.olfEnabled
        if state.phys.mcAcq.olfShuntEnabled
            state.phys.mcAcq.olfValveStatus = 0;
            mcAcqValveSwitch(state.phys.mcAcq.olfDefaultValve, 0)
            stop(state.phys.mcAcq.olfDevice);            
        else
            state.phys.mcAcq.olfValveStatus = 0;            
            mcAcqValveSwitch(state.phys.mcAcq.olfDefaultValve)
            stop(state.phys.mcAcq.olfDevice);
        end
    end
    % END MOD