function resetSpines
global gh state

for handle=state.imageProc.spine.linehandles
	try
		delete(handle)
	end
end

for handle=state.imageProc.spine.texthandles
	try
		delete(handle)
	end
end

state.imageProc.spine.linehandles=[];
state.imageProc.spine.texthandles=[];
state.imageProc.spine.spineLengths=[];


state.imageProc.spine.numberOfSpines = 0;
updateGUIByGlobal('state.imageProc.spine.numberOfSpines');

state.imageProc.spine.spineDensity = 0;
updateGUIByGlobal('state.imageProc.spine.spineDensity');
