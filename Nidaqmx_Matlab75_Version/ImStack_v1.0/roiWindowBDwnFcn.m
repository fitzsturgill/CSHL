function roiWindowBDwnFcn
% Function for moving a patch object interactively...
global state gh
type = get(gh.roiAnalysis.figure1,'SelectionType');
switch type
case 'extend'	%shift click
	set(gh.roiAnalysis.figure1, 'WindowButtonMotionFcn', 'movePatchObject');
end
