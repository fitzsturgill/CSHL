function updateFullFileName(h)
	global state
	if ~ischar(state.files.baseName)
		state.files.baseName='';
	end
	if ~ischar(state.files.savePath)
		state.files.savePath='';
	end
	try
		state.files.fullFileName=fullfile(state.files.savePath, [state.files.baseName zeroPadNum2str(state.files.fileCounter)]);
	catch
	end
	