function readMax(fname, pname)
	global state
	state.internal.lastTaskDone=3;
	state.internal.status=0;
	try
		imageInfo=imfinfo(fullfile(pname, [fname '.tif']));
		length(imageInfo);
		state.lastHeader=imageInfo(1).ImageDescription;
	catch
		lasterr
		return
	end
	
	pixelsPerLine=state.acq.pixelsPerLine
	valueFromHeaderString(state.acq.maxData=cell(1,3);
	try
		for counter=1:state.init.maximumNumberOfInputChannels
			if valueFromHeaderString(['state.acq.savingChannel' num2str(counter)], state.lastHeader)
				counter
				state.acq.maxData{counter}=imread(fullfile(pname, [fname '.tif']), counter);
				set(state.internal.maximagehandle(counter), 'EraseMode', 'none', 'CData', state.acq.maxData{counter});
			else
				state.acq.maxData{counter}=[];
				set(state.internal.maximagehandle(counter), 'EraseMode', 'none', 'CData', state.acq.maxData{counter});
			end
		end
	catch
		lasterr
	end
			
	