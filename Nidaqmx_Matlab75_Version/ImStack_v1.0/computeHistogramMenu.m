function computeHistogramMenu
global state gh

	state.imageProc.internal.histFigure = figure('NumberTitle', 'off', ...
		'Name', 'Histogram','Color','White');
	value = get(gh.imageProcessingGUI.fileName, 'Value');
	[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});

	colormap(gray);
	cl = class(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.cell.currentFrame{value}));
	
	switch cl
	case 'uint8'	
		imhist(double(state.imageProc.cell.currentImage{value}(y:y1,x:x1,state.imageProc.cell.currentFrame{value}))/double(state.imageProc.highPixelValue));
	case 'uint16'
		imhist(double(state.imageProc.cell.currentImage{value}(y:y1,x:x1,state.imageProc.cell.currentFrame{value}))/double(state.imageProc.highPixelValue));
	otherwise
		imhist(state.imageProc.cell.currentImage{value}(y:y1,x:x1,state.imageProc.cell.currentFrame{value}));
	end
