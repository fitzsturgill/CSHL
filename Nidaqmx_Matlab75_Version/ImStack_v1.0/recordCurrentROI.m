function recordCurrentROI
% This function will add the current axis limits to the list of ROI's
global state
[x,x1,y,y1]=getCurrentAxisLimits(gca);
state.imageProc.ROIPositionVector=[y y1 x x1];
disp(['Recorded Position ' num2str(state.imageProc.ROIPositionVector) ' as ROI.']);
