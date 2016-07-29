function figureButtonOverCallbackImagePr
global state gh

% This function will grab the current position and intensity from the figure
% and display it in the imageGUI Window.

	currentPoint = recordCurrentPoint(gca);
	state.imageProc.internal.positionX = currentPoint(1,1);
	updateGUIByGlobal('state.imageProc.internal.positionX');
	state.imageProc.internal.positionY = currentPoint(1,2);
	updateGUIByGlobal('state.imageProc.internal.positionY');
	CData = get(gco, 'CData');
	state.imageProc.internal.intensity = CData(state.imageProc.internal.positionY, state.imageProc.internal.positionX);
	updateGUIByGlobal('state.imageProc.internal.intensity');
    [x,x1,y,y1] = getCurrentAxisLimits(get(gco,'Parent'));
    state.imageProc.internal.meanIntensity = mean2(CData(y:y1,x:x1));
	updateGUIByGlobal('state.imageProc.internal.meanIntensity');
	state.imageProc.internal.sumIntensity = sum(sum((CData(y:y1,x:x1))));
	updateGUIByGlobal('state.imageProc.internal.sumIntensity');
	figureButtonDownFcn;