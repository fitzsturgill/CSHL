function savePulseSet
% writes pulseSet to disk

	global state
	
	if isempty(state.phys.pulses.pulseSetPath) | isempty(state.phys.pulses.pulseSetName)
		savePulseSetAs;
	else
		fn=fieldnames(state.phys.pulses);
		
		pulseSet=[];
		
		for counter=1:length(fn)
			if findstr('List', fn{counter})
				eval(['pulseSet.' fn{counter} '=state.phys.pulses. ' fn{counter} ';']);
			end
		end

		save(fullfile(state.phys.pulses.pulseSetPath, state.phys.pulses.pulseSetName), 'pulseSet')
		%save(fullfile(state.phys.pulses.pulseSetPath, state.phys.pulses.pulseSetName), 'pulseSet', '-v6')

        setPhysStatusString('pulseSet saved');
        state.phys.pulses.pulseSetChanged=0;
	end
