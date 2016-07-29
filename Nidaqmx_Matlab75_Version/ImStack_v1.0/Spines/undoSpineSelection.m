function undoSpineSelection
	global state
	if state.imageProc.spine.numberOfSpines<=1
		state.imageProc.spine.numberOfSpines = 0;
		state.imageProc.spine.spineLengths = [];
	else
		state.imageProc.spine.numberOfSpines = state.imageProc.spine.numberOfSpines - 1;
		state.imageProc.spine.spineLengths = state.imageProc.spine.spineLengths(1:state.imageProc.spine.numberOfSpines);
	end
		updateGUIByGlobal('state.imageProc.spine.numberOfSpines');
	end
	saveSpineDataToExcel;