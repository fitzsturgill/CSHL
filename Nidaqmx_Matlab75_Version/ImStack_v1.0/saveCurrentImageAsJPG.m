function saveCurrentImageAsJPG
global state gh
	
	try
		cd(state.imageProc.internal.savePath);
	catch
	end
	
	[fname, pname] = uiputfile({'*.jpg'},'Save Region of Interest as...' );
	
	if strcmp(state.imageProc.internal.savePath, 'C:\matlabR12\work')
		
		state.imageProc.internal.savePath = pname;
		updateGUIByGlobal('state.imageProc.internal.savePath');
		state.imageProc.internal.saveBaseName = fname;
		updateGUIByGlobal('state.imageProc.internal.saveBaseName');
		
	elseif ~strcmp(state.imageProc.internal.savePath, 'C:\matlabR12\work') & ...
			state.imageProc.internal.updateSavePath == 1
	
		state.imageProc.internal.savePath = pname;
		updateGUIByGlobal('state.imageProc.internal.savePath');
		state.imageProc.internal.saveBaseName = fname;
		updateGUIByGlobal('state.imageProc.internal.saveBaseName');
	
	else
	end
	
	[path, name, ext] = fileparts(state.imageProc.fileName);
	cd(path);
	value = get(gh.imageProcessingGUI.fileName, 'Value');
	[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
	% convert picture to uint 8 
	jpgarray = state.imageProc.cell.currentImage{value}(y:y1, x:x1,state.imageProc.cell.currentFrame{value});
	maxpixel = max(max(jpgarray));
	jpgarray = double(jpgarray)*(255/double(maxpixel));
	jpgarray = uint8(jpgarray);
	if ischar(fname)		
		imwrite(jpgarray, [pname fname '.jpg']);	
	else
		return
	end








