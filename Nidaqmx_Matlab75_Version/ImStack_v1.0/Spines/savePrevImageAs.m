function savePrevImageAs
global gh state

[fname, pname] = uiputfile('*.tif','Save Preview Image As...' );
	if fname < 1
		return
	else
		cd(pname);
		printHandleToFile(gh.spineGUI.previewaxis, fname);
	end
