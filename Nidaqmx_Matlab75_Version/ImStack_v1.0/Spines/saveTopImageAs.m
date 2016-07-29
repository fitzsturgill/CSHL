function saveTopImageAs
global gh state

[fname, pname] = uiputfile('*.tif','Save Top Image As...' );
	if fname < 1
		return
	else
		cd(pname);
		printHandleToFile(gh.spineGUI.initialaxis, fname);
	end
	