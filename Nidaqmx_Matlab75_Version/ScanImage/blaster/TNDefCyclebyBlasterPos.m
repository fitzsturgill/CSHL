function TNDefCyclebyBlasterPos()
	global state
  
    state.cycle.imageOnList(:) = 0;
    state.cycle.physOnList(:) = 0;
    
    for pos=1:length(state.blaster.indexList)
          
		cyc=2*pos-1; %track image trials

		state.cycle.imageOnList(cyc) = 1;
		state.cycle.repeatsList(cyc) = 1;
		state.cycle.delayList(cyc) = 3;
		state.cycle.framesList(cyc) = 1;
		state.cycle.numberOfZSlicesList(cyc) = 1;
		state.cycle.zStepSizeList(cyc) = 1;
		state.cycle.blasterList(cyc) = 0;
		state.cycle.trackerList(cyc) = 1;
		state.cycle.avgFramesList(cyc) = 0;
		state.cycle.linescanList(cyc) = 0;
		state.cycle.physOnList(cyc) = 0;
		state.cycle.stagePosList(cyc) = 0;
		state.cycle.scanSetupList(cyc) = 0;
		state.cycle.da0List(cyc) = 0;
		state.cycle.da1List(cyc) = 0;
		state.cycle.aux4List(cyc) = 0;
		state.cycle.aux5List(cyc) = 0;
		state.cycle.aux6List(cyc) = 0;
		state.cycle.aux7List(cyc) = 0;
			
		cyc=2*pos;%line scan trials

		state.cycle.imageOnList(cyc) = 1;
		state.cycle.repeatsList(cyc) = 1;
		state.cycle.delayList(cyc) = 5;
		state.cycle.framesList(cyc) = 1;
		state.cycle.numberOfZSlicesList(cyc) = 1;
		state.cycle.zStepSizeList(cyc) = 1;
		state.cycle.blasterList(cyc) = pos;
		state.cycle.trackerList(cyc) = 0;
		state.cycle.avgFramesList(cyc) = 0;
		state.cycle.linescanList(cyc) = 0;
		state.cycle.physOnList(cyc) = 1;
		state.cycle.stagePosList(cyc) = 0;
		state.cycle.scanSetupList(cyc) = 0;
		state.cycle.da0List(cyc) = 19;
		state.cycle.da1List(cyc) = 0;
		state.cycle.aux4List(cyc) = 0;
		state.cycle.aux5List(cyc) = 0;
		state.cycle.aux6List(cyc) = 0;
		state.cycle.aux7List(cyc) = 0;

		
    end
	
	state.internal.cycleChanged=1;
	updateCycleDisplay(1);
