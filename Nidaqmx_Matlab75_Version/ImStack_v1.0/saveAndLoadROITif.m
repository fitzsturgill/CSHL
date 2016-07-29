function saveAndLoadROITif
global state gh

	[path,name,ext] = fileparts(state.imageProc.fileName);
	try
		cd(path);
	catch
	end
	
	[fname, pname] = uiputfile('*.tif','Save Region of Interest as...' );

	value = get(gh.imageProcessingGUI.fileName, 'Value');
	header = get(gh.imageParsingGUI.header, 'String')
	[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
	
	if fname > 0
			arrayToTiff(state.imageProc.cell.currentImage{value} ...
				(y:y1, x:x1,state.imageProc.cell.currentFrame{value}:state.imageProc.cell.numberOfFrames{value}), [pname fname], header);
			loadImageFromName([pname fname '.tif'], path);	
	else
		return
	end
	


