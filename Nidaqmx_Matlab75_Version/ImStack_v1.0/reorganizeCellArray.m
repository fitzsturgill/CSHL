function cellOutput = reorganizeCellArrays
global state gh

oldString = get(gh.imageProcessingGUI.fileName, 'String');
oldValue = get(gh.imageProcessingGUI.fileName, 'Value');

if ischar(oldString)
	set(gh.imageProcessingGUI.fileName, 'String', 'Load Image...');

		removeFieldFromCellArray(state.imageProc.cell.currentImage, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.fileName, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.numberOfFrames, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.pixelsPerLine, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.linesPerFrame, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.montageStart, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.montageEnd, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.numberOfChannels, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.numberOfFrames, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.scanRotation, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.averaged, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.numberOfZSlices, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.maxEnd, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.maxStart, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.urrentFrame, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.lowPixelValue, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.highPixelValue, oldValue,1);
		
		close(state.imageProc.internal.Figure);
		removeFieldFromCellArray(state.imageProc.internal.Figure, oldValue,1);
		removeFieldFromCellArray(state.imageProc.internal.axis, oldValue,1);
		removeFieldFromCellArray(state.imageProc.internal.imagehandle, oldValue,1);
		
else
	newstring = removeFieldFromCellArray(oldstring, oldValue, 0);
	set(gh.imageProcessingGUI.fileName, 'String', 'Load Image...');
	
	removeFieldFromCellArray(state.imageProc.cell.currentImage, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.fileName, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.numberOfFrames, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.pixelsPerLine, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.linesPerFrame, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.montageStart, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.montageEnd, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.numberOfChannels, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.numberOfFrames, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.scanRotation, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.averaged, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.parsing.numberOfZSlices, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.maxEnd, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.maxStart, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.urrentFrame, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.lowPixelValue, oldValue,1);
		removeFieldFromCellArray(state.imageProc.cell.highPixelValue, oldValue,1);
		
		close(state.imageProc.internal.Figure);
		removeFieldFromCellArray(state.imageProc.internal.Figure, oldValue,1);
		removeFieldFromCellArray(state.imageProc.internal.axis, oldValue,1);
		removeFieldFromCellArray(state.imageProc.internal.imagehandle, oldValue,1);		
end
state.imageProc.internal.imageCounter = state.imageProc.internal.imageCounter - 1;
	
	