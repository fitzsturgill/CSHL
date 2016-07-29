function newBox
	global state
	
	state.pcell.currentStartX=0;
	state.pcell.currentStartY=0;
	state.pcell.currentEndX=0;
	state.pcell.currentEndY=0;
	state.pcell.currentActiveStatus=1;
	state.pcell.currentPowerLevel=0;
	state.pcell.currentFrameNumber=1;
	
	updateGuiByGlobal('state.pcell.currentStartX');
	updateGuiByGlobal('state.pcell.currentStartY');
	updateGuiByGlobal('state.pcell.currentEndX');
	updateGuiByGlobal('state.pcell.currentEndY');
	updateGuiByGlobal('state.pcell.currentActiveStatus');
	updateGuiByGlobal('state.pcell.currentPowerLevel');
	updateGuiByGlobal('state.pcell.currentFrameNumber');

	if ishandle(state.pcell.currentBoxHandle)
		delete(state.pcell.currentBoxHandle);
		state.pcell.currentBoxHandle=-1;
	end
	storeBoxSettings;

	