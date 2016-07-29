function incrementExcelRow
	global state
	state.imageProc.spine.row = state.imageProc.spine.row+1;
	updateGUIByGlobal('state.imageProc.spine.row');

