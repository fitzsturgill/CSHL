function AutoSelectFiles
global gh state

try 
	cd(state.imageProc.pathName);
catch
	state.imageProc.pathName = pwd;
end
state.imageProc.autoSelect = 1;
updateGUIByGlobal('state.imageProc.autoSelect');
state.imageProc.names = selectFilesFromList(state.imageProc.pathName);
