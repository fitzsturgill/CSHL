function loadROIIntoSpine
global state gh

	value = get(gh.imageProcessingGUI.fileName, 'Value');
	[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
	eval(['	state.imageProc.spine.initialImage  '...
			' = state.imageProc.cell.currentImage{value}(y:y1, x:x1,state.imageProc.cell.currentFrame{value}:state.imageProc.cell.numberOfFrames{value});']);
	
	columns = size(state.imageProc.spine.initialImage,2);
	rows = size(state.imageProc.spine.initialImage,1);
	
	state.imageProc.spine.lowPixelValue = min(min(state.imageProc.spine.initialImage));
	updateGUIByGlobal('state.imageProc.spine.lowPixelValue');
	state.imageProc.spine.highPixelValue = .4*double(max(max(state.imageProc.spine.initialImage)));
	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
	
	set(gh.spineGUI.figure1,'Colormap', makeColormap('gray',8));
	set(gh.spineGUI.initialaxis,'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue], 'Parent', gh.spineGUI.figure1);
	state.imageProc.spine.imagehandle = image('CData', state.imageProc.spine.initialImage, ...
		'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis, 'ButtonDownFcn', 'spineImageOverFcn');
	seeGUI('gh.spineGUI.figure1');