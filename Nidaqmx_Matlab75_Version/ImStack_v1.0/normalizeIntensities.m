function normalizeIntensities(image)
global gh state

	state.imageProc.internal.histEqFigure = figure('NumberTitle', 'off', ...
		'Name', 'Equalized Image');
	value = get(gh.imageProcessingGUI.fileName, 'Value');
	colormap(gray);
	histeq(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.cell.currentFrame{value}));
