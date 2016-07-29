function clearSpineStats
global gh state

state.imageProc.spine.spineCounter = 1;
state.imageProc.spine.globalSpineArea =0;
state.imageProc.spine.globalSpineLength=0;
state.imageProc.spine.totalSpines = 0;
updateGUIByGlobal('state.imageProc.spine.totalSpines');
state.imageProc.spine.totaldendriteLength = 0;
updateGUIByGlobal('state.imageProc.spine.totaldendriteLength');
state.imageProc.spine.totalDensity = 0;
updateGUIByGlobal('state.imageProc.spine.totalDensity');