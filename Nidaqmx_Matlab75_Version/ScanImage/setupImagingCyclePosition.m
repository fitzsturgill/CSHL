function setupImagingCyclePosition
	global state

	%disp(['Calling setupimagingCyclePosition' ]);
	
	redo=0;
	needToAllocate=0;
	
	state.acq.zStepSize = state.cycle.zStepSizeList(state.cycle.currentCyclePosition);
	updateGuiByGlobal('state.acq.zStepSize');
	
	if state.acq.autoTrack ~= state.cycle.trackerList(state.cycle.currentCyclePosition)
		state.acq.autoTrack=state.cycle.trackerList(state.cycle.currentCyclePosition);
		updateGuiByGlobal('state.acq.autoTrack');
	end
		
	if state.cycle.blasterList(state.cycle.currentCyclePosition)>0
		if ~state.blaster.active
			state.blaster.active=1;
			updateGuiByGlobal('state.blaster.active');
			if state.blaster.blankImaging
				state.internal.needNewPcellPowerOutput=1;
			end	
			redo=1;
		end
		if state.blaster.currentConfig ~= state.cycle.blasterList(state.cycle.currentCyclePosition)
			redo=1;
			global gh
			c=get(gh.blaster.Config, 'Children');
			set(c, 'Checked', 'off');
			f=find(strcmp(get(c, 'Label'), [num2str(state.cycle.blasterList(state.cycle.currentCyclePosition))...
					' ' state.blaster.allConfigs{state.cycle.blasterList(state.cycle.currentCyclePosition),1}]));
			if ~isempty(f)
				set(c(f(1)), 'Checked', 'on');
			end
		end
	elseif state.blaster.active
		state.blaster.active=0;
		updateGuiByGlobal('state.blaster.active');
		if state.blaster.blankImaging
			state.internal.needNewPcellPowerOutput=1;
		end	
		redo=1;
	end
	state.blaster.currentConfig = state.cycle.blasterList(state.cycle.currentCyclePosition);
		
	%% Added by TN Apr 01 2005
	try
		if state.cycle.scanSetupList(state.cycle.currentCyclePosition) ~= 0 
			restoreScan(state.cycle.scanSetupList(state.cycle.currentCyclePosition));
		end
	catch
		disp(['ERROR ' ]);
	end
	%% End TN Apr 01 2005
	
	if state.acq.lineScanCB ~= state.cycle.linescanList(state.cycle.currentCyclePosition)
		state.acq.lineScanCB = state.cycle.linescanList(state.cycle.currentCyclePosition);
		state.internal.needNewRotatedMirrorOutput=1;
		updateGuiByGlobal('state.acq.lineScanCB');
		redo=1;
	end
	
	if state.acq.numberOfZSlices ~= state.cycle.numberOfZSlicesList(state.cycle.currentCyclePosition)
		state.acq.numberOfZSlices = state.cycle.numberOfZSlicesList(state.cycle.currentCyclePosition);
		updateGuiByGlobal('state.acq.numberOfZSlices');
		needToAllocate=1;
	end
	
	if state.acq.numberOfFrames ~= state.cycle.framesList(state.cycle.currentCyclePosition)
		state.acq.numberOfFrames = state.cycle.framesList(state.cycle.currentCyclePosition);
		updateGuiByGlobal('state.acq.numberOfFrames');
		set(state.daq.grabInput, 'SamplesPerTrigger', state.internal.samplesPerFrame*state.acq.numberOfFrames);
		needToAllocate=1;
		redo=1;
	end
	
	if state.acq.averaging~=state.cycle.avgFramesList(state.cycle.currentCyclePosition)
		state.acq.averaging=state.cycle.avgFramesList(state.cycle.currentCyclePosition);
		updateGuiByGlobal('state.acq.averaging');
		needToAllocate=1;
	end

	if needToAllocate
		preallocateMemory;
	end

	if redo
		state.internal.needNewRepeatedMirrorOutput=1;
		state.internal.needNewPcellRepeatedOutput=1;
	end
	


	
	