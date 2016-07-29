function timerUserSettings_Physiology
	global state

	loadPulseSet(state.phys.pulses.pulseSetPath, state.phys.pulses.pulseSetName);

	changeChannelType(0);
	changeChannelType(1);
	
	if ~state.analysisMode
		state.phys.settings.outputRate=setverify(state.phys.daq.outputDevice, 'SampleRate', state.phys.settings.outputRate*1000)/1000;
		state.phys.settings.inputRate=setverify(state.phys.daq.inputDevice, 'SampleRate', state.phys.settings.inputRate*1000)/1000;;
	end
	updateGuiByGlobal('state.phys.settings.outputRate');
	updateGuiByGlobal('state.phys.settings.inputRate');

	setupPhysDaqInputChannels;
