function chooseLoadPath
global gh state

try 
	cd(state.imageProc.pathName);
catch
end

[name, path] = uigetfile({'*.tif;*.cfd;*.CFD'}, 'Choose a File in the Load path and hit "Open"...');
if name == 0
	return
else
	set(gh.autotransformGUI.pathLoadName, 'String', path);
	genericCallback(gh.autotransformGUI.pathLoadName);
	fullname = [path name];
	set(gh.autotransformGUI.baseLoadName, 'String', getbasename(fullname,path));
	genericCallback(gh.autotransformGUI.baseLoadName);
end

state.imageProc.acquisitionNumber = 1;
updateGUIByGlobal('state.imageProc.acquisitionNumber');


		