function make3DHistogramMenu
global state gh

	value = get(gh.imageProcessingGUI.fileName, 'Value');

	state.imageProc.internal.hist3DFigure = figure('NumberTitle', 'off', ...
		'Name', '3D Contour');
	
	value = get(gh.imageProcessingGUI.fileName, 'Value');
	[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
	
	surfc(double(state.imageProc.cell.currentImage{value}(y:y1,x:x1,state.imageProc.cell.currentFrame{value})));
