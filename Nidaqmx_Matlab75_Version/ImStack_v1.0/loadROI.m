function loadROI
global state gh

	value = get(gh.imageProcessingGUI.fileName, 'Value');
	[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
	eval(['	state.imageProc.internal.ROI' num2str(state.imageProc.internal.ROICounter) ...
			' = state.imageProc.cell.currentImage{value}(y:y1, x:x1,state.imageProc.cell.currentFrame{value}:state.imageProc.cell.numberOfFrames{value});']);
	loadImageFromArray(['state.imageProc.internal.ROI' num2str(state.imageProc.internal.ROICounter)]);
	state.imageProc.internal.ROICounter = state.imageProc.internal.ROICounter+1;
	

