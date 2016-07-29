function chooseSavePath
global gh state

try
	cd(state.imageProc.internal.savePath);
catch
end

[name, path] = uiputfile([state.imageProc.baseName '00.tif'], 'Choose the Save path and hit "SAVE"...');
if name == 0
	return
else
	set(gh.autotransformGUI.pathSaveName, 'String', path);
	genericCallback(gh.autotransformGUI.pathSaveName);
	fullname = [path name];
	set(gh.autotransformGUI.baseSaveName, 'String', getbasename(fullname,path));
	genericCallback(gh.autotransformGUI.baseSaveName);
end
