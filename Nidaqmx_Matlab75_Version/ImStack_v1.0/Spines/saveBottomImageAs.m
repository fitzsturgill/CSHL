function saveBottomImageAs
global gh state

[fname, pname] = uiputfile('*.tif','Save Bottom Image As...' );
	if fname < 1
		return
	else
		cd(pname);
		printHandleToFile(gh.spineGUI.initialaxis2, fname);
	end
