function timerStart_Physiology
	global state gh
	timerSetPackageStatus(1, 'Physiology');

	state.cycle.lastPulseUsed0 = state.cycle.pulseToUse0;
	state.cycle.lastPulseUsed1 = state.cycle.pulseToUse1;
    if state.phys.daq.auxOutputBoardIndex
    	state.cycle.lastUsedAuxPulses = [state.cycle.aux4List(state.cycle.currentCyclePosition) ...
                state.cycle.aux5List(state.cycle.currentCyclePosition) ...
                state.cycle.aux6List(state.cycle.currentCyclePosition) ...
                state.cycle.aux7List(state.cycle.currentCyclePosition)];
    end
	
	state.files.lastAcquisition=state.files.fileCounter;

	set(gh.physControls.startButton, 'String', 'ABORT');
	set(gh.scope.start, 'Enable', 'off');
  
  
    %if(state.cycle.imageOnList(state.cycle.currentCyclePosition))   % TN 09May05
    %    start([state.phys.daq.inputDevice state.phys.daq.outputDevice state.phys.daq.auxOutputDevice]);
    %else
        start([state.phys.daq.inputDevice state.phys.daq.outputDevice]);
    %end
    
	setPhysStatusString('Running...');
	if state.phys.daq.auxOutputBoardIndex & any(state.cycle.lastUsedAuxPulses) & ~state.cycle.imageOnList(state.cycle.currentCyclePosition)
		start(state.phys.daq.auxOutputDevice);
    end	
    
    %FS MOD
    if ~isempty(state.phys.mcAcq.mcInputBoardIndex) && state.cycle.mcOnList(state.cycle.currentCyclePosition)
        start(state.phys.mcAcq.mcInputDevice);
    end
    %
