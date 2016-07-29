function setupScanningPcellDAQ
	global state
	
	if state.pcell.pcellOn == 1
		% for grab %
		state.init.ao1 = analogoutput('nidaq',state.pcell.pcellBoardIndex);			 % 1 refers to the NIDAQ MIO64 Board
		for pcellCounter=1:state.pcell.numberOfPcells
			eval(['state.init.pockellChannel' num2str(pcellCounter) ...
					'= addchannel(state.init.ao1, state.pcell.pcellChannelIndex' num2str(pcellCounter) ');']);	
		end
		set(state.init.ao1, 'SampleRate', state.acq.outputRate);
		set(state.init.ao1, 'TriggerType', 'HWDigital');	
		syncNIDAQBoards(state.init.ao2, state.init.ao1);	% Will sync the 2 board clocks.

		% for focus %
		state.init.ao1F = analogoutput('nidaq',state.init.pockellBoardIndex);			 % 1 refers to the NIDAQ MIO64 Board
		for pcellCounter=1:state.pcell.numberOfPcells
			eval(['state.init.pockellChannel' num2str(pcellCounter) ...
					'= addchannel(state.init.ao1F, state.pcell.pcellChannelIndex' num2str(pcellCounter) ');']);	
		end
		set(state.init.ao1F, 'SampleRate', state.acq.outputRate);
		set(state.init.ao1F, 'TriggerType', 'HWDigital');	
	end
