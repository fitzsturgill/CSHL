function getROIFromImage(x,x1,y,y1)
%Loads ROI from current frame by indices given.
% Use this to get limis after click on the image : [x,x1,y,y1]=getCurrentAxisLimits(gca)
global state
state.imageProc.ROISelected=state.imageProc.currentImage(y:y1,x:x1,state.imageProc.currentFrame:state.imageProc.numberOfFrames);
loadImageFromArray('state.imageProc.ROISelected');