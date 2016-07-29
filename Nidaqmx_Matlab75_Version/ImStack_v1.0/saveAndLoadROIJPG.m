function saveAndLoadROIJPG
global state gh

	[path,name,ext] = fileparts(state.imageProc.fileName);
	try
		cd(path);
	catch
	end
	
	[fname, pname] = uiputfile('*.jpg','Save Region of Interest as...' );

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
	
	if ischar(fname)
		fullname=[pname fname '.jpg'];
		state.imageProc.currentjpeg = [];
		state.imageProc.currentjpeg = imread(fullname);
		size(state.imageProc.currentjpeg)
		state.imageProc.currentjpeg = rgb2gray(state.imageProc.currentjpeg);
		loadImageFromArray('state.imageProc.currentjpeg');
		set(gh.fileCounterGUI.filetype, 'Value', 2);
		cd(pname);
	else
		return
	end



