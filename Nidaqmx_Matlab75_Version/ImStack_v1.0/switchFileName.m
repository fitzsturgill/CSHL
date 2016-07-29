function switchFileName(Value);
global state gh

% This is the callback to select the current image characteristics to be displayed.

		% set the first current image array to the first image opened.
		updateGUIByGlobalCell('state.imageProc.currentImage', Value);
		
		state.imageProc.fileNameGUI = Value;
		updateGUIByGlobal('state.imageProc.fileNameGUI');
		
		% update header
		updateGUIByGlobalCell('state.imageProc.parsing.header', Value);
		
		% get the size of the image (total number of frames)
		updateGUIByGlobalCell('state.imageProc.numberOfFrames', Value);
		
		% get the pixels per line (X) of the image 
		updateGUIByGlobalCell('state.imageProc.parsing.pixelsPerLine', Value);
		
		% get the lines per frame (Y)of the image 
		updateGUIByGlobalCell('state.imageProc.parsing.linesPerFrame', Value);
	
		% PArsing Updates
		updateGUIByGlobalCell('state.imageProc.parsing.numberOfChannels', Value);
		updateGUIByGlobalCell('state.imageProc.parsing.numberOfFrames', Value);
		updateGUIByGlobalCell('state.imageProc.parsing.scanRotation', Value);
		updateGUIByGlobalCell('state.imageProc.parsing.averaged', Value);
		updateGUIByGlobalCell('state.imageProc.parsing.numberOfZSlices', Value);
		
		% update montage end
		updateGUIByGlobalCell('state.imageProc.montageStart', Value);
		updateGUIByGlobalCell('state.imageProc.montageEnd', Value);
		
		% update max end
		updateGUIByGlobalCell('state.imageProc.maxEnd', Value);
		updateGUIByGlobalCell('state.imageProc.maxStart', Value);
		
		% update movie end
		updateGUIByGlobalCell('state.imageProc.movieEnd', Value);
		updateGUIByGlobalCell('state.imageProc.movieStart', Value);
		
		% update the current frame in the image
		updateGUIByGlobalCell('state.imageProc.currentFrame', Value);
	
		% Set lookup table to scale for image
		updateGUIByGlobalCell('state.imageProc.lowPixelValue', Value);
		updateGUIByGlobalCell('state.imageProc.highPixelValue', Value);
	
		% Make this figure Current if it exists
		try
			set(state.imageProc.internal.Figure{Value}, 'Visible', 'On')
		catch
		end

		updateSliderMaxMin;