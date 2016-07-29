function saveAndLoadROI
global state gh

	[path,name,ext] = fileparts(state.imageProc.fileName);
	try
		cd(path);
	catch
	end
	
	[fname, pname] = uiputfile('*.tif','Save Region of Interest as...' );

	value = get(gh.imageProcessingGUI.fileName, 'Value');
	[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});

	if ischar(fname)
		arrayToTiff(state.imageProc.cell.currentImage{value} ...
			(y:y1, x:x1,state.imageProc.cell.currentFrame{value}:state.imageProc.cell.numberOfFrames{value}), [pname fname], '');
		loadImageFromName([pname fname '.tif'], path);
	
	else
		return
	end



