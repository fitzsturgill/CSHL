function inputCurrentROI
% This function will add the entered values to the list of ROI's
global state

prompt  = {'Start Y:','End Y:','Start X:','End X:'};
title   = 'Manual Selection of ROI';
lines= 1;
if ~isempty(state.imageProc.ROIPositionVector)
def     = {num2str(state.imageProc.ROIPositionVector(1)),num2str(state.imageProc.ROIPositionVector(2)),...
	num2str(state.imageProc.ROIPositionVector(3)),num2str(state.imageProc.ROIPositionVector(4))};
else
	def={'1','2','3','4'};
end

answer  = inputdlg(prompt,title,lines,def);
if ~isempty(answer)
	state.imageProc.ROIPositionVector(1)=str2num(answer{1});
	state.imageProc.ROIPositionVector(2)=str2num(answer{2});
	state.imageProc.ROIPositionVector(3)=str2num(answer{3});
	state.imageProc.ROIPositionVector(4)=str2num(answer{4});
	disp(['Recorded Position ' num2str(state.imageProc.ROIPositionVector) ' as ROI.']);

end
