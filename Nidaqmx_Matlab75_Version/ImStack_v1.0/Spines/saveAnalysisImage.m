function saveAnalysisImage
	global state gh
	
	if state.imageProc.spine.topImage
		name = state.imageProc.spine.loadedFileNameTop;
	elseif state.imageProc.spine.bottomImage
		name = state.imageProc.spine.loadedFileNameBot;
	else
		name =''
	end
	
	saveas(figure(get(gh.spineGUI.mainAxes, 'parent')), fullfile(state.imageProc.spine.openPath, [name '_spines.fig']));
	