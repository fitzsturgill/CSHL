function saveMaxProjAxis(handle)
global gh state

try 
	cd(state.imageProc.internal.savePath);
end

[fname, pname] = uiputfile('*.tif', 'Save File As...');

if fname < 1
	return
end

filename = [ pname fname];
printHandleToFile(handle, filename);
state.imageProc.internal.savePath = pname;
updateGUIByGlobal('state.imageProc.internal.savePath');