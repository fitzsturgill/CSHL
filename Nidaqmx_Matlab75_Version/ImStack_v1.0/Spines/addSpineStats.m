function addSpineStats
global gh state

if state.imageProc.spine.spineCounter == 1
	state.imageProc.spine.globalSpineArea = [state.imageProc.spine.spineFeatures.Area];
	state.imageProc.spine.globalSpineLength = [state.imageProc.spine.spineFeatures.MajorAxisLength];
	state.imageProc.spine.totalSpines = length(state.imageProc.spine.globalSpineArea);
	updateGUIByGlobal('state.imageProc.spine.totalSpines');
	state.imageProc.spine.totaldendriteLength = state.imageProc.spine.dendriteLength;
	updateGUIByGlobal('state.imageProc.spine.totaldendriteLength');
	state.imageProc.spine.totalDensity = state.imageProc.spine.totalSpines/state.imageProc.spine.dendriteLength;
	updateGUIByGlobal('state.imageProc.spine.totalDensity');
	
	state.imageProc.spine.spineCounter = state.imageProc.spine.spineCounter +1;
else
	state.imageProc.spine.globalSpineArea = [state.imageProc.spine.globalSpineArea state.imageProc.spine.spineFeatures.Area];
	state.imageProc.spine.globalSpineLength = [state.imageProc.spine.globalSpineLength state.imageProc.spine.spineFeatures.MajorAxisLength];
	state.imageProc.spine.totalSpines = length(state.imageProc.spine.globalSpineArea);
	updateGUIByGlobal('state.imageProc.spine.totalSpines');
	state.imageProc.spine.totaldendriteLength = state.imageProc.spine.totaldendriteLength + state.imageProc.spine.dendriteLength;
	updateGUIByGlobal('state.imageProc.spine.totaldendriteLength');
	state.imageProc.spine.totalDensity = state.imageProc.spine.totalSpines/state.imageProc.spine.totaldendriteLength;
	updateGUIByGlobal('state.imageProc.spine.totalDensity');
	state.imageProc.spine.spineCounter = state.imageProc.spine.spineCounter +1;
end
