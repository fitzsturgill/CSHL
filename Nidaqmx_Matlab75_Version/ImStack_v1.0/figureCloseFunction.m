function figureCloseFunction 
global state gh

try
	
name = get(gcf, 'Name');
oldValue = getMenuIndex2(gh.imageProcessingGUI.fileName, name);

if ischar(oldValue)
	oldValue = 1;
end

oldString = get(gh.imageProcessingGUI.fileName, 'String');

if oldValue == 1
	newValue = oldValue;
else
	newValue = oldValue-1;
end

if ischar(oldString) 
	set(gh.imageProcessingGUI.fileName, 'Value', newValue);
	set(gh.maxProjectionGUI.fileName, 'Value', newValue);
	set(gh.montageGUI.fileName, 'Value', newValue);
	set(gh.movieGUI.fileName, 'Value', newValue);
	set(gh.overlayGUI.fileName1, 'Value', newValue); 
	set(gh.overlayGUI.fileName2, 'Value', newValue); 
	set(gh.overlayGUI.fileName3, 'Value', newValue); 
	set(gh.averagingGUI.fileName, 'Value', newValue); 
	set(gh.mathGUI.fileName1, 'Value', newValue); 
	set(gh.mathGUI.fileName2, 'Value', newValue); 
	set(gh.roiAnalysis.roiPopupmenu, 'Value', newValue);
    
	set(gh.imageProcessingGUI.fileName, 'String', 'Load Image...');
	set(gh.montageGUI.fileName, 'String', 'Load Image...');
	set(gh.movieGUI.fileName, 'String', 'Load Image...');
	set(gh.maxProjectionGUI.fileName, 'String', 'Load Image...');
	set(gh.overlayGUI.fileName1, 'String', 'Load Image...'); 
	set(gh.overlayGUI.fileName2, 'String', 'Load Image...'); 
	set(gh.overlayGUI.fileName3, 'String', 'Load Image...'); 
	set(gh.averagingGUI.fileName, 'String', 'Load Image...');
	set(gh.mathGUI.fileName1, 'String', 'Load Image...'); 
	set(gh.mathGUI.fileName2, 'String', 'Load Image...'); 
	set(gh.imageParsingGUI.filename, 'String', 'Load Image...');
    set(gh.roiAnalysis.roiPopupmenu, 'String', 'Load Image...');
	
		state.imageProc.parsing.cell.header= removeFieldFromCellArray(state.imageProc.parsing.cell.header, oldValue,1);
		state.imageProc.cell.currentImage = removeFieldFromCellArray(state.imageProc.cell.currentImage, oldValue,1);
		state.imageProc.cell.fileName = removeFieldFromCellArray(state.imageProc.cell.fileName, oldValue,1);
		state.imageProc.cell.numberOfFrames = removeFieldFromCellArray(state.imageProc.cell.numberOfFrames, oldValue,1);
		state.imageProc.parsing.cell.pixelsPerLine = removeFieldFromCellArray(state.imageProc.parsing.cell.pixelsPerLine, oldValue,1);
		state.imageProc.parsing.cell.linesPerFrame = removeFieldFromCellArray(state.imageProc.parsing.cell.linesPerFrame, oldValue,1);
		state.imageProc.cell.montageStart = removeFieldFromCellArray(state.imageProc.cell.montageStart, oldValue,1);
		state.imageProc.cell.montageEnd = removeFieldFromCellArray(state.imageProc.cell.montageEnd, oldValue,1);
		state.imageProc.parsing.cell.numberOfChannels = removeFieldFromCellArray(state.imageProc.parsing.cell.numberOfChannels, oldValue,1);
		state.imageProc.parsing.cell.numberOfFrames = removeFieldFromCellArray(state.imageProc.parsing.cell.numberOfFrames, oldValue,1);
		state.imageProc.parsing.cell.scanRotation = removeFieldFromCellArray(state.imageProc.parsing.cell.scanRotation, oldValue,1);
		state.imageProc.parsing.cell.averaged = removeFieldFromCellArray(state.imageProc.parsing.cell.averaged, oldValue,1);
		state.imageProc.parsing.cell.numberOfZSlices = removeFieldFromCellArray(state.imageProc.parsing.cell.numberOfZSlices, oldValue,1);
		state.imageProc.cell.maxEnd = removeFieldFromCellArray(state.imageProc.cell.maxEnd, oldValue,1);
		state.imageProc.cell.maxStart = removeFieldFromCellArray(state.imageProc.cell.maxStart, oldValue,1);
		state.imageProc.cell.movieEnd = removeFieldFromCellArray(state.imageProc.cell.movieEnd, oldValue,1);
		state.imageProc.cell.movieStart = removeFieldFromCellArray(state.imageProc.cell.movieStart, oldValue,1);
		state.imageProc.cell.currentFrame = removeFieldFromCellArray(state.imageProc.cell.currentFrame, oldValue,1);
		state.imageProc.cell.lowPixelValue = removeFieldFromCellArray(state.imageProc.cell.lowPixelValue, oldValue,1);
		state.imageProc.cell.highPixelValue = removeFieldFromCellArray(state.imageProc.cell.highPixelValue, oldValue,1);
		
		set(state.imageProc.internal.Figure{oldValue}, 'Visible', 'off');
		state.imageProc.internal.Figure = removeFieldFromCellArray(state.imageProc.internal.Figure, oldValue,1);
		state.imageProc.internal.axis = removeFieldFromCellArray(state.imageProc.internal.axis, oldValue,1);
		state.imageProc.internal.imagehandle = removeFieldFromCellArray(state.imageProc.internal.imagehandle, oldValue,1);
		
		resetFieldsToDefault;
		
elseif length(oldString) == 1
		set(gh.imageProcessingGUI.fileName, 'Value', newValue);
		set(gh.maxProjectionGUI.fileName, 'Value', newValue);
		set(gh.montageGUI.fileName, 'Value', newValue);
		set(gh.movieGUI.fileName, 'Value', newValue);
		set(gh.overlayGUI.fileName1, 'Value', newValue); 
		set(gh.overlayGUI.fileName2, 'Value', newValue); 
		set(gh.overlayGUI.fileName3, 'Value', newValue); 
		set(gh.averagingGUI.fileName, 'Value', newValue); 
		set(gh.mathGUI.fileName1, 'Value', newValue); 
		set(gh.mathGUI.fileName2, 'Value', newValue); 
	    set(gh.roiAnalysis.roiPopupmenu, 'Value', newValue);
        
		set(gh.imageProcessingGUI.fileName, 'String', 'Load Image...');
		set(gh.montageGUI.fileName, 'String', 'Load Image...');
		set(gh.movieGUI.fileName, 'String', 'Load Image...');
		set(gh.maxProjectionGUI.fileName, 'String', 'Load Image...');
		set(gh.overlayGUI.fileName1, 'String', 'Load Image...'); 
		set(gh.overlayGUI.fileName2, 'String', 'Load Image...'); 
		set(gh.overlayGUI.fileName3, 'String', 'Load Image...'); 
		set(gh.averagingGUI.fileName, 'String', 'Load Image...');			
		set(gh.imageParsingGUI.filename, 'String', 'Load Image...');
		set(gh.mathGUI.fileName1, 'String', 'Load Image...'); 
		set(gh.mathGUI.fileName2, 'String', 'Load Image...'); 
        set(gh.roiAnalysis.roiPopupmenu, 'String', 'Load Image...');
         
		state.imageProc.parsing.cell.header= removeFieldFromCellArray(state.imageProc.parsing.cell.header, oldValue,1);
		state.imageProc.cell.currentImage = removeFieldFromCellArray(state.imageProc.cell.currentImage, oldValue,1);
		state.imageProc.cell.fileName = removeFieldFromCellArray(state.imageProc.cell.fileName, oldValue,1);
		state.imageProc.cell.numberOfFrames = removeFieldFromCellArray(state.imageProc.cell.numberOfFrames, oldValue,1);
		state.imageProc.parsing.cell.pixelsPerLine = removeFieldFromCellArray(state.imageProc.parsing.cell.pixelsPerLine, oldValue,1);
		state.imageProc.parsing.cell.linesPerFrame = removeFieldFromCellArray(state.imageProc.parsing.cell.linesPerFrame, oldValue,1);
		state.imageProc.cell.montageStart = removeFieldFromCellArray(state.imageProc.cell.montageStart, oldValue,1);
		state.imageProc.cell.montageEnd = removeFieldFromCellArray(state.imageProc.cell.montageEnd, oldValue,1);
		state.imageProc.cell.movieStart = removeFieldFromCellArray(state.imageProc.cell.movieStart, oldValue,1);
		state.imageProc.cell.movieEnd = removeFieldFromCellArray(state.imageProc.cell.movieEnd, oldValue,1);
		state.imageProc.parsing.cell.numberOfChannels = removeFieldFromCellArray(state.imageProc.parsing.cell.numberOfChannels, oldValue,1);
		state.imageProc.parsing.cell.numberOfFrames = removeFieldFromCellArray(state.imageProc.parsing.cell.numberOfFrames, oldValue,1);
		state.imageProc.parsing.cell.scanRotation = removeFieldFromCellArray(state.imageProc.parsing.cell.scanRotation, oldValue,1);
		state.imageProc.parsing.cell.averaged = removeFieldFromCellArray(state.imageProc.parsing.cell.averaged, oldValue,1);
		state.imageProc.parsing.cell.numberOfZSlices = removeFieldFromCellArray(state.imageProc.parsing.cell.numberOfZSlices, oldValue,1);
		state.imageProc.cell.maxEnd = removeFieldFromCellArray(state.imageProc.cell.maxEnd, oldValue,1);
		state.imageProc.cell.maxStart = removeFieldFromCellArray(state.imageProc.cell.maxStart, oldValue,1);
		state.imageProc.cell.currentFrame = removeFieldFromCellArray(state.imageProc.cell.currentFrame, oldValue,1);
		state.imageProc.cell.lowPixelValue = removeFieldFromCellArray(state.imageProc.cell.lowPixelValue, oldValue,1);
		state.imageProc.cell.highPixelValue = removeFieldFromCellArray(state.imageProc.cell.highPixelValue, oldValue,1);
		
		set(state.imageProc.internal.Figure{oldValue}, 'Visible', 'off');
		state.imageProc.internal.Figure = removeFieldFromCellArray(state.imageProc.internal.Figure, oldValue,1);
		state.imageProc.internal.axis = removeFieldFromCellArray(state.imageProc.internal.axis, oldValue,1);
		state.imageProc.internal.imagehandle = removeFieldFromCellArray(state.imageProc.internal.imagehandle, oldValue,1);
		
		resetFieldsToDefault;
		
elseif length(oldString) > 1
	newstring = removeFieldFromCellArray(oldString, oldValue, 0);
	
	set(gh.imageProcessingGUI.fileName, 'String', newstring);
	set(gh.montageGUI.fileName, 'String', newstring);
	set(gh.maxProjectionGUI.fileName, 'String', newstring);
	set(gh.movieGUI.fileName, 'String', newstring);
	set(gh.overlayGUI.fileName1, 'String', newstring); 
	set(gh.overlayGUI.fileName2, 'String', newstring); 
	set(gh.overlayGUI.fileName3, 'String', newstring); 
	set(gh.averagingGUI.fileName, 'String', newstring);
	set(gh.imageParsingGUI.filename, 'String', newstring{newValue});
	set(gh.mathGUI.fileName1, 'String', newstring); 
	set(gh.mathGUI.fileName2, 'String', newstring); 
    set(gh.roiAnalysis.roiPopupmenu, 'String', newstring);
	
	set(gh.imageProcessingGUI.fileName, 'Value', newValue);
	set(gh.maxProjectionGUI.fileName, 'Value', newValue);
	set(gh.montageGUI.fileName, 'Value', newValue);
	set(gh.movieGUI.fileName, 'Value', newValue);
	set(gh.overlayGUI.fileName1, 'Value', newValue); 
	set(gh.overlayGUI.fileName2, 'Value', newValue); 
	set(gh.overlayGUI.fileName3, 'Value', newValue); 
	set(gh.averagingGUI.fileName, 'Value', newValue);
	set(gh.mathGUI.fileName1, 'Value', newValue); 
	set(gh.mathGUI.fileName2, 'Value', newValue); 
     set(gh.roiAnalysis.roiPopupmenu, 'Value', newValue);
    
	
		state.imageProc.parsing.cell.header= removeFieldFromCellArray(state.imageProc.parsing.cell.header, oldValue,1);
		state.imageProc.cell.currentImage = removeFieldFromCellArray(state.imageProc.cell.currentImage, oldValue,1);
		state.imageProc.cell.fileName = removeFieldFromCellArray(state.imageProc.cell.fileName, oldValue,1);
		state.imageProc.cell.numberOfFrames = removeFieldFromCellArray(state.imageProc.cell.numberOfFrames, oldValue,1);
		state.imageProc.parsing.cell.pixelsPerLine = removeFieldFromCellArray(state.imageProc.parsing.cell.pixelsPerLine, oldValue,1);
		state.imageProc.parsing.cell.linesPerFrame = removeFieldFromCellArray(state.imageProc.parsing.cell.linesPerFrame, oldValue,1);
		state.imageProc.cell.montageStart = removeFieldFromCellArray(state.imageProc.cell.montageStart, oldValue,1);
		state.imageProc.cell.montageEnd = removeFieldFromCellArray(state.imageProc.cell.montageEnd, oldValue,1);
		state.imageProc.cell.movieStart = removeFieldFromCellArray(state.imageProc.cell.movieStart, oldValue,1);
		state.imageProc.cell.movieEnd = removeFieldFromCellArray(state.imageProc.cell.movieEnd, oldValue,1);
		state.imageProc.parsing.cell.numberOfChannels = removeFieldFromCellArray(state.imageProc.parsing.cell.numberOfChannels, oldValue,1);
		state.imageProc.parsing.cell.numberOfFrames = removeFieldFromCellArray(state.imageProc.parsing.cell.numberOfFrames, oldValue,1);
		state.imageProc.parsing.cell.scanRotation = removeFieldFromCellArray(state.imageProc.parsing.cell.scanRotation, oldValue,1);
		state.imageProc.parsing.cell.averaged = removeFieldFromCellArray(state.imageProc.parsing.cell.averaged, oldValue,1);
		state.imageProc.parsing.cell.numberOfZSlices = removeFieldFromCellArray(state.imageProc.parsing.cell.numberOfZSlices, oldValue,1);
		state.imageProc.cell.maxEnd = removeFieldFromCellArray(state.imageProc.cell.maxEnd, oldValue,1);
		state.imageProc.cell.maxStart = removeFieldFromCellArray(state.imageProc.cell.maxStart, oldValue,1);
		state.imageProc.cell.currentFrame = removeFieldFromCellArray(state.imageProc.cell.currentFrame, oldValue,1);
		state.imageProc.cell.lowPixelValue = removeFieldFromCellArray(state.imageProc.cell.lowPixelValue, oldValue,1);
		state.imageProc.cell.highPixelValue = removeFieldFromCellArray(state.imageProc.cell.highPixelValue, oldValue,1);
		
		set(state.imageProc.internal.Figure{oldValue}, 'Visible', 'off');
		state.imageProc.internal.Figure = removeFieldFromCellArray(state.imageProc.internal.Figure, oldValue,1);
		state.imageProc.internal.axis = removeFieldFromCellArray(state.imageProc.internal.axis, oldValue,1);
		state.imageProc.internal.imagehandle = removeFieldFromCellArray(state.imageProc.internal.imagehandle, oldValue,1);
		switchFileName(newValue);
		state.imageProc.internal.imageCounter = state.imageProc.internal.imageCounter - 1;

end

catch
end

	
	
	