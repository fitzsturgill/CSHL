function movePatchObject
% This will move the patch object specified to the current axis location.
global gh state
currentPoint=get(gh.roiAnalysis.roiAxes,'CurrentPoint');
currentPoint=currentPoint(1,1:2);

if ishandle(state.imageProc.roiAnalysis.patchHandle)
	%Follow current Point with the mouse...
	state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
	CurrentroiXpos=mean(state.imageProc.roiAnalysis.XPatch);
	difference=currentPoint(1)-CurrentroiXpos;
	set(state.imageProc.roiAnalysis.patchHandle,'XData',state.imageProc.roiAnalysis.XPatch+difference);
	state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
	state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
	CurrentroiYpos=mean(state.imageProc.roiAnalysis.YPatch);
	difference=currentPoint(2)-CurrentroiYpos;
	set(state.imageProc.roiAnalysis.patchHandle,'YData',state.imageProc.roiAnalysis.YPatch+difference);
	state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
	% Reset mover positions....
	state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
	state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
	state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
	state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
	state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
	updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
	updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');
end

