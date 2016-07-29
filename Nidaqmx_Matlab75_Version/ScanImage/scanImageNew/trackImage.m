function trackImage
	global state gh

	if ~state.acq.autoTrack 
		return
	end

	if ~state.internal.keepAllSlicesInMemory & (state.acq.numberOfZSlices>1)
		beep;
		disp('*** Stack based image track skipped ***');
		disp('*** Requires keepAllSlicesInMemory to be activated from the settings menu ***');
		return
	end

	try
		setStatusString('Finding shifts ...');

		state.acq.shifts=zeros(3, state.acq.numberOfZSlices);
		maxCC=-1000;
		for sliceCounter=1:state.acq.numberOfZSlices
			if state.acq.averaging
				startSlice=sliceCounter;
				stopSlice=sliceCounter;
			else
				startSlice=(sliceCounter-1)*state.acq.numberOfFrames+1;
				stopSlice=startSlice+state.acq.numberOfFrames-1;
			end

			state.acq.trackerImage=mean(state.acq.acquiredData{state.acq.trackerChannel}(:,:,startSlice:stopSlice),3);
			state.acq.shifts(:, sliceCounter)...
				=findShift(state.acq.trackerReference, medfilt2(state.acq.acquiredData{state.acq.trackerChannel}(:,:,sliceCounter)));
            disp(['CC ' num2str(state.acq.shifts(3, sliceCounter))])
            if state.acq.shifts(3, sliceCounter)>maxCC
				maxCC=state.acq.shifts(3, sliceCounter);
				maxCCSlice=sliceCounter;
				shift=state.acq.shifts(1:2, sliceCounter);
			end
        end
        
        addEntryToNotebook(2, ['trackImage acq ' num2str(state.files.fileCounter) ' (' ... 
            num2str(shift(2)) ', ' num2str(shift(1)) ', ' num2str(maxCCSlice) ') CC=' num2str(maxCC)]);



		state.acq.pixelShiftY=shift(1)-state.internal.trackerY0+1;
		state.acq.pixelShiftX=shift(2)-state.internal.trackerX0+1;
		state.acq.sliceShift=state.acq.zStepSize*(maxCCSlice-1);
		
		updateGUIByGlobal('state.acq.pixelShiftX');
		updateGUIByGlobal('state.acq.pixelShiftY');
		updateGUIByGlobal('state.acq.sliceShift');
		scanShiftX=2*state.acq.scanAmplitudeX*(state.acq.pixelShiftX/state.acq.pixelsPerLine)/state.acq.zoomFactor;
		scanShiftY=2*state.acq.scanAmplitudeY*(state.acq.pixelShiftY/state.acq.linesPerFrame)/state.acq.zoomFactor;

		c = cos(state.acq.scanRotation*pi/180);
		s = sin(state.acq.scanRotation*pi/180);
		state.acq.scanShiftX =  c * scanShiftX + s * scanShiftY;
		state.acq.scanShiftY = -s * scanShiftX + c * scanShiftY;				

		if (state.acq.pixelShiftY > state.acq.maxShift) | (state.acq.pixelShiftX > state.acq.maxShift)
			beep;
			setStatusString('SHIFT TOO LARGE');
			disp('*** UNABLE to apply shift. Shift too large.');
            state.piezo.next_pos=state.internal.saveScanInfo(state.cycle.scanSetupList(state.cycle.currentCyclePosition), 13);
            piezoUpdatePositionNoWait;
		else
			state.acq.postRotOffsetX=state.acq.postRotOffsetX+state.acq.scanShiftX;
			state.acq.postRotOffsetY=state.acq.postRotOffsetY+state.acq.scanShiftY;		

			if ~isempty(state.blaster.XList)
				state.blaster.XList=state.blaster.XList+state.acq.scanShiftX;
				state.blaster.YList=state.blaster.YList+state.acq.scanShiftY;
				state.blaster.X=state.blaster.XList(state.blaster.displayPos);
				state.blaster.Y=state.blaster.YList(state.blaster.displayPos);
				updateGuiByGlobal('state.blaster.X');
				updateGuiByGlobal('state.blaster.Y');
			end
			state.internal.refShiftX = state.internal.refShiftX + state.acq.scanShiftX;
			state.internal.refShiftY = state.internal.refShiftY + state.acq.scanShiftY;
			state.internal.needNewRotatedMirrorOutput=1;
			state.internal.needNewPcellRepeatedOutput=1;
			
			%TN 01Apr05
			if(state.cycle.scanSetupList(state.cycle.currentCyclePosition)>0)
				updateSavedScans(shift, state.cycle.scanSetupList(state.cycle.currentCyclePosition));
				if state.acq.numberOfZSlices==3
					switch maxCCSlice
						%update savedscaninfo and move there
						case 1
							state.internal.saveScanInfo(state.cycle.scanSetupList(state.cycle.currentCyclePosition), 13)=state.internal.saveScanInfo(state.cycle.scanSetupList(state.cycle.currentCyclePosition), 13)-state.cycle.zStepSize;
                            disp('Moving Up');
                        case 3
                            state.internal.saveScanInfo(state.cycle.scanSetupList(state.cycle.currentCyclePosition), 13)=state.internal.saveScanInfo(state.cycle.scanSetupList(state.cycle.currentCyclePosition), 13)+state.cycle.zStepSize;						
                            disp('Moving Down');
                    end
                    %restoreScan(state.cycle.scanSetup)
					state.piezo.next_pos=state.internal.saveScanInfo(state.cycle.scanSetupList(state.cycle.currentCyclePosition), 13);
					piezoUpdatePositionNoWait;
				end
 			end
		end
 	catch
 		setStatusString('AUTOTRACK ERROR');
 		disp('*** trackImage: Error in autoTrack  ');
 		disp(lasterr);
 	end

	