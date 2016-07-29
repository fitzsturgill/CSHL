function varargout=zoomCurrentROI
% This function will zoom current axis to saved ROI the current saved axis limits

out=0;
global state
if ~isempty(state.imageProc.ROIPositionVector)
   set(gca,'XLim',[state.imageProc.ROIPositionVector(3) state.imageProc.ROIPositionVector(4)],...
	   'YLim',[state.imageProc.ROIPositionVector(1) state.imageProc.ROIPositionVector(2)])
   out=1;
end

if nargout==1
    varargout{1}=out;
end
