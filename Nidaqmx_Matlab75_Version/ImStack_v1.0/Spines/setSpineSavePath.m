function setSpineSavePath
global gh state

[name, path] = uiputfile('Choose A Path.tif', 'Choose the Save path and hit "SAVE"...');
if isnumeric(path)
	return
end

state.imageProc.spine.pathName = path;
updateGUIByGlobal('state.imageProc.spine.pathName');
