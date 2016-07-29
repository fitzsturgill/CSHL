function roiPatchButtonDwn
% This is the button down functiuon for the patch object used in the ROI Analysis
global gh state
type = get(gh.roiAnalysis.figure1,'SelectionType');
switch type
case 'open'	%Double Click
case 'normal'	% Do the analysis
	state.imageProc.roiAnalysis.roiAngle=state.imageProc.roiAnalysis.roiAngle-10;
	updateGUIByGlobal('state.imageProc.roiAnalysis.roiAngle');
	global gh state
	% Rotate object about its center....
	if ishandle(state.imageProc.roiAnalysis.patchHandle)
		state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
		state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
		state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
		state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
		state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
		updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
		updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');
		
		rotate(state.imageProc.roiAnalysis.patchHandle,[0 0 1],...
        state.imageProc.roiAnalysis.roiAngle-state.imageProc.roiAnalysis.lastAngle,state.imageProc.roiAnalysis.patchCenter );

		state.imageProc.roiAnalysis.lastAngle=state.imageProc.roiAnalysis.roiAngle;
		
		state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
		state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
		state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
		state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
		state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
		updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
		updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');
	end
	roiStats;
	
	
case 'alt'	% Right click
	state.imageProc.roiAnalysis.roiAngle=state.imageProc.roiAnalysis.roiAngle+10;
	updateGUIByGlobal('state.imageProc.roiAnalysis.roiAngle');
	global gh state
	% Rotate object about its center....
	if ishandle(state.imageProc.roiAnalysis.patchHandle)
		state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
		state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
		state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
		state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
		state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
		updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
		updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');
		
		rotate(state.imageProc.roiAnalysis.patchHandle,[0 0 1],...
        state.imageProc.roiAnalysis.roiAngle-state.imageProc.roiAnalysis.lastAngle,state.imageProc.roiAnalysis.patchCenter );

		state.imageProc.roiAnalysis.lastAngle=state.imageProc.roiAnalysis.roiAngle;
		
		state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
		state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
		state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
		state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
		state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
		updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
		updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');
	end
		roiStats;
case 'extend'
end