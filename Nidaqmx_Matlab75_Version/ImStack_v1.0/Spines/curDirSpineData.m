function curDirSpineData
global gh state

[name, path] = uigetfile({'*.txt'}, 'Choose a File in the Load path and hit "Open"...');
if name == 0
	return
else
	state.imageProc.spineData.openPath = path;
	cd(path);
end

