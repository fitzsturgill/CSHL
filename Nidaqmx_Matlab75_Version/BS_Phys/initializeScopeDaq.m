function initializeScopeDaq

	global state
	if state.analysisMode
		return
	end
	
	state.phys.daq.scopeOutputDevice = analogoutput('nidaq', state.phys.daq.outputBoardIndex);
	set(state.phys.daq.scopeOutputDevice, 'SampleRate', state.phys.scope.outputRate*1000);
	state.phys.scope.actualOutputRate=get(state.phys.daq.scopeOutputDevice, 'SampleRate')/1000;
	set(state.phys.daq.scopeOutputDevice, 'TriggerType', 'HwDigital');		

	state.phys.daq.scopeInputDevice = analoginput('nidaq', state.phys.daq.inputBoardIndex);
	set(state.phys.daq.scopeInputDevice, 'TriggerType', 'HwDigital');				
	set(state.phys.daq.scopeInputDevice, 'SampleRate', state.phys.scope.inputRate*1000);
    set(state.phys.daq.scopeInputDevice, 'InputType', 'SingleEnded'); %FS MOD
	state.phys.scope.actualInputRate=get(state.phys.daq.scopeInputDevice, 'SampleRate')/1000;
	set(state.phys.daq.scopeInputDevice, 'SamplesAcquiredFcn', {'processScopeData'});
