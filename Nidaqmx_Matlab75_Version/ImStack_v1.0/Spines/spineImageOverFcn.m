function spineImageOverFcn
global gh state

if state.imageProc.spine.topImage
		currentPoint = recordCurrentPoint(gh.spineGUI.initialaxis);
		state.imageProc.spine.xpos = currentPoint(1,1);
		updateGUIByGlobal('state.imageProc.spine.xpos');
		state.imageProc.spine.ypos = currentPoint(1,2);
		updateGUIByGlobal('state.imageProc.spine.ypos');
	if state.imageProc.spine.maxFlag == 0
		state.imageProc.spine.intensity = state.imageProc.spine.initialImage(state.imageProc.spine.ypos, state.imageProc.spine.xpos, ...
			state.imageProc.spine.currentSpineFrame);
		updateGUIByGlobal('state.imageProc.spine.intensity');
	else
		state.imageProc.spine.intensity = state.imageProc.spine.maxProjection(state.imageProc.spine.ypos, state.imageProc.spine.xpos);
		updateGUIByGlobal('state.imageProc.spine.intensity');
	end

else
	return
end
