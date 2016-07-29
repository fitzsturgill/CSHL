function loadPreviewSpine
global gh state

[fname, pname] = uigetfile('*.tif', 'Choose image to load as Preview');

	if ~isnumeric(fname)
		periods=findstr(fname, '.');
		if any(periods)								
			fname=fname(1:periods(1)-1);
		else
			disp('Found File without Extension');
			return
		end	
		
		% Compute the full file name
		filename = fullfile(pname, [fname '.tif']);
		cd(pname);
	else
		return
	end
	
	state.imageProc.spine.previewImageData = opentif(filename);
	rows = size(state.imageProc.spine.previewImageData,1);
	columns = size(state.imageProc.spine.previewImageData,2);
	state.imageProc.spine.XLim3 = [1 columns];
	state.imageProc.spine.YLim3 = [1 rows];
	set(gh.spineGUI.previewaxis, 'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'Clim', [state.imageProc.spine.lowPixelValue ...
			state.imageProc.spine.highPixelValue]);
	try
		state.imageProc.spine.previewImage = image('CData', state.imageProc.spine.previewImageData,...
		'CDataMapping', 'scaled','Parent', gh.spineGUI.previewaxis);
	catch
		state.imageProc.spine.previewImage = image('CData', state.imageProc.spine.previewImageData(:,:,1),...
			'CDataMapping', 'scaled','Parent', gh.spineGUI.previewaxis);
		disp('File Multi-Tiff; Only Showing First Frame');
	end

	set(gh.spineGUI.previewFigure, 'visible', 'on');