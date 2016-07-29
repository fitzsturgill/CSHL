function spineImageOverFcnBot
global gh state

if state.imageProc.spine.bottomImage
		currentPoint = recordCurrentPoint(gh.spineGUI.initialaxis2);
		state.imageProc.spine.xpos = currentPoint(1,1);
		updateGUIByGlobal('state.imageProc.spine.xpos');
		state.imageProc.spine.ypos = currentPoint(1,2);
		updateGUIByGlobal('state.imageProc.spine.ypos');
		
		if state.imageProc.spine.maxFlag2 == 0
			
			state.imageProc.spine.intensity = state.imageProc.spine.initialImage2(state.imageProc.spine.ypos, state.imageProc.spine.xpos, ...
				state.imageProc.spine.currentSpineFrame2);
			updateGUIByGlobal('state.imageProc.spine.intensity');
		else
			state.imageProc.spine.intensity = state.imageProc.spine.maxProjection2(state.imageProc.spine.ypos, state.imageProc.spine.xpos);
			updateGUIByGlobal('state.imageProc.spine.intensity');
		end
		
else
	return
end
