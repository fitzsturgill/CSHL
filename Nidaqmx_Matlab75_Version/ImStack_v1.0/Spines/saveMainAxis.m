function saveMainAxis
global gh state

if state.imageProc.spine.autoSave == 0
	
	[fname, pname] = uiputfile('*.tif','Save Analyzed Image as...' );
	if fname < 1
		return
	else
		cd(pname);
		printHandleToFile(gh.spineGUI.mainAxes, fname);
	end
	if state.imageProc.spine.savePreview
			[fname, pname] = uiputfile('*.tif','Save Preview Image as...' );
			if isnumeric(pname)
				return
			end
			printHandleToFile(gh.spineGUI.previewaxis, fname);
	end
else
	if strcmp(state.imageProc.spine.pathName, 'C:\matlabR12\work\')
		button = questdlg('Would You Like to Choose a Save Path Now?','No Save Path Selected',...
		'Yes','Use Default','Cancel', 'Yes');
		if strcmp(button,'Yes')
 		  setSpineSavePath;
		  cd(state.imageProc.spine.pathName);
	  	elseif strcmp(button,'Use Default')
		  disp('Using Default Path');
		  cd(state.imageProc.spine.pathName);
		else
			return
		end
	else
		cd(state.imageProc.spine.pathName);
	end
	
	if state.imageProc.spine.fileSave < 10
		fname = [state.imageProc.spine.baseName '0' num2str(state.imageProc.spine.fileSave) '.tif'];
	else
		fname = [state.imageProc.spine.baseName   num2str(state.imageProc.spine.fileSave) '.tif'];
	end
		printHandleToFile(gh.spineGUI.mainAxes, fname);
		if state.imageProc.spine.savePreview
			
			if state.imageProc.spine.fileSave < 10
				fname = [state.imageProc.spine.baseName(1:(end-1)) 'p00' num2str(state.imageProc.spine.fileSave) '.tif'];
			else
				fname = [state.imageProc.spine.baseName(1:(end-1)) 'p0' num2str(state.imageProc.spine.fileSave) '.tif'];
			end
			printHandleToFile(gh.spineGUI.previewaxis, fname);
		end
		state.imageProc.spine.fileSave = state.imageProc.spine.fileSave+1;
		updateGUIByGlobal('state.imageProc.spine.fileSave');
end

		