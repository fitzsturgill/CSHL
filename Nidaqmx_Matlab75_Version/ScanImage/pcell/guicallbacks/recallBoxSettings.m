function recallBoxSettings(boxNum)
	global state
	
	if nargin<1
		boxNum=state.pcell.currentBoxNumber;
	end
	
	state.pcell.currentStartX=state.pcell.boxListStartX(boxNum);
	state.pcell.currentStartY=state.pcell.boxListStartY(boxNum);
	state.pcell.currentEndX=state.pcell.boxListEndX(boxNum);
	state.pcell.currentEndY=state.pcell.boxListEndY(boxNum);
	state.pcell.currentActiveStatus=state.pcell.boxListActive(boxNum);
	state.pcell.currentBoxHandle=state.pcell.boxListHandles(boxNum);
	state.pcell.currentPowerLevel=state.pcell.boxListPowerLevel(state.pcell.currentBoxNumber);
	state.pcell.currentBoxChannel=state.pcell.boxListBoxChannel(state.pcell.currentBoxNumber);
	state.pcell.currentFrameNumber=state.pcell.boxListFrameNumber(state.pcell.currentBoxNumber);

% 	updateGuiByGlobal('state.pcell.currentStartX');
% 	updateGuiByGlobal('state.pcell.currentStartY');
% 	updateGuiByGlobal('state.pcell.currentEndX');
% 	updateGuiByGlobal('state.pcell.currentEndY');
	updateGuiByGlobal('state.pcell.currentActiveStatus');
	updateGuiByGlobal('state.pcell.currentPowerLevel');
	updateGuiByGlobal('state.pcell.currentBoxChannel');
	updateGuiByGlobal('state.pcell.currentFrameNumber');

	if ishandle(state.pcell.currentBoxHandle) & state.pcell.currentBoxHandle~=0
		set(state.pcell.currentBoxHandle,'Color',[0 0 1]);
		drawnow;
	end

	
