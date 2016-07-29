function resetDendrite
global gh state 

for handle=state.imageProc.spine.dendriteLines
	try
		 delete(handle);
	catch
	end
end
state.imageProc.spine.dendriteLines = [];

state.imageProc.spine.dendriteLength = 0;
updateGUIByGlobal('state.imageProc.spine.dendriteLength');

state.imageProc.spine.spineDensity = 0;
updateGUIByGlobal('state.imageProc.spine.spineDensity');