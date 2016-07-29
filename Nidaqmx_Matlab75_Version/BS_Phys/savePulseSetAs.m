function savePulseSetAs
% save pulse set to disk with new name

	global state
	
	if ~isempty(state.phys.pulses.pulseSetPath)
		try
			cd(state.phys.pulses.pulseSetPath)
		catch
		end
	end
	
	[fname, pname]=uiputfile('*.mat', 'Choose file for pulse set');

	if ~isnumeric(fname)
		periods=findstr(fname, '.');
		if any(periods)								
			fname=fname(1:periods(1)-1);
		end		
		state.phys.pulses.pulseSetName = fname;
		state.phys.pulses.pulseSetPath = pname;
		updateGUIByGlobal('state.phys.pulses.pulseSetName');
		savePulseSet;
	else
		setPhysStatusString('Cannot open file');
	end
	