function endAcquisition
	global state gh

	% endAcquisition.m*****
	% Function called at the end of the acquistion that will park the laser, close the shutter,
	% write the data to disk, reset the counters (internal), reset the currentMode, and make the 
	% Grab One and Loop buttons visible.
	%
	% Written By: Thomas Pologruto and Bernardo Sabatini
	% Cold Spring Harbor Labs
	% March 2, 2001
	
	% Code for displaying Max Projections
	setPcellsToDefault;
	if (state.acq.numberOfFrames == 1 | state.acq.averaging == 1) & any(state.acq.maxImage)
		if state.internal.keepAllSlicesInMemory % BSMOD 1/18/2
			position = state.internal.zSliceCounter + 1;
		else
			position = 1;
		end

        if state.acq.dualLaserMode==1   % lasers go simulataneously therefor one image window per laser
            channelList=1:state.init.maximumNumberOfInputChannels;
        else
            channelList=[1:state.init.maximumNumberOfInputChannels 11:10+state.init.maximumNumberOfInputChannels];
        end

		for channelCounter = channelList
	        inputChannel=mod(channelCounter, 10);		
			if getfield(state.acq, ['maxImage' num2str(inputChannel)]) ...		% If max is on and ...
				& getfield(state.acq,['acquiringChannel' num2str(inputChannel)])	% channel is on
				if state.internal.zSliceCounter==0			% BSMOD2 2/27/2
					if state.acq.maxMode==0
						state.acq.maxData{channelCounter} = state.acq.acquiredData{channelCounter}(:,:,position);
					else
						state.acq.maxData{channelCounter} = state.acq.acquiredData{channelCounter}(:,:,position);
					end
				else
					if state.acq.maxMode==0
						state.acq.maxData{channelCounter} = max(state.acq.acquiredData{channelCounter}(:,:,position), ...
							state.acq.maxData{channelCounter});
					else
						state.acq.maxData{channelCounter} = ...
							(state.acq.acquiredData{channelCounter}(:,:,state.internal.zSliceCounter + 1) + ...
							state.internal.zSliceCounter*state.acq.maxData{channelCounter})/(state.internal.zSliceCounter + 1);	
						%  BSMOD 1/18/2 eliminated reliance on position for above 2 lines
					end					
				end
					% Displays the current Max images on the screen as they are acquired.
				set(state.internal.maximagehandle(channelCounter), 'EraseMode', 'none', 'CData', ...
					state.acq.maxData{channelCounter}); 	
			end
				
		end
		drawnow;	
	end
	
	if state.internal.zSliceCounter + 1 == state.acq.numberOfZSlices
	% Done Acquisition.
		if state.files.autoSave		% BSMOD - Check status of autoSave option
			setStatusString('Writing data...');
			writeData;
			writeMaxData;
		end
		
		parkMirrors;
		
		state.internal.zSliceCounter = state.internal.zSliceCounter + 1;
		updateGUIByGlobal('state.internal.zSliceCounter');
			
		
		if state.acq.numberOfZSlices > 1
			if state.piezo.usePiezo
				wait(state.piezo.Output, state.piezo.tsec+1) ;
			else
				mp285FinishMove(1);	% check that movement worked during stack
			end
			executeGoHome;
		end				

		timerRegisterPackageDone('Imaging');	
	elseif state.internal.zSliceCounter < state.acq.numberOfZSlices - 1
	% Between Acquisitions or ZSlices
		setStatusString('Next Slice...');
		
		if state.files.autoSave		% BSMOD - Check status of autoSave option
			setStatusString('Writing data...');
			writeData;
		end
	
		state.internal.zSliceCounter = state.internal.zSliceCounter + 1;
		updateGUIByGlobal('state.internal.zSliceCounter');
	
		state.internal.frameCounter = 0;
		updateGUIByGlobal('state.internal.frameCounter');
				
		% if there is slice specific pcell control, remake the pcell output
		if state.pcell.boxSliceSpecific
			makeNewPcellRepeatedOutput;
		end

		setStatusString('Acquiring...');
	
		putDataGrab;

		if state.piezo.usePiezo
			wait(state.piezo.Output, state.piezo.tsec+1) ;
		else
			mp285FinishMove(0);	% check that movement worked during stack
		end
		
		startGrab;
		openShutter;
		diotrigger;
	end
	

	
	
	