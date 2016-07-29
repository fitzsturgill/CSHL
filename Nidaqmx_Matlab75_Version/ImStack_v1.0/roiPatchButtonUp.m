function roiPatchButtonUp
% This is the button down functiuon for the patch object used in the ROI Analysis
global gh state
type = get(gh.roiAnalysis.figure1,'SelectionType');
switch type
case 'open'	%Double Click
case 'normal'	% Do the analysis
	set(gh.roiAnalysis.figure1, 'WindowButtonMotionFcn', '');	
case 'alt'	% Right click
	
case 'extend'
end