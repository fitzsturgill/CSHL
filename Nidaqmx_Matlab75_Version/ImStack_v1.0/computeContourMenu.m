function computeContourMenu
global state gh


	state.imageProc.internal.contourFigure = figure('NumberTitle', 'off', ...
		'Name', 'Contour');
	value = get(gh.imageProcessingGUI.fileName, 'Value');
	[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
	
	colormap(gray);
	imcontour(state.imageProc.cell.currentImage{value}(y:y1,x:x1,state.imageProc.cell.currentFrame{value}));
	