function varargout=loadCurrentROI
% This function will load the current saved axis limits

% % Make UI Context Menus
% cmenuHandle=createUiContextMenu('Save ROI', 'recordCurrentROI', 'Load Saved ROI','loadCurrentROI','Load and Process ROI',...
%     'loadAndProcessCurROI','Median Filter', 'applyMedianFilter');
% addUimenuToHandle(state.imageProc.internal.imagehandle{state.imageProc.internal.imageCounter}, cmenuHandle);
% 
out=0;
global state
if ~isempty(state.imageProc.ROIPositionVector)
    state.imageProc.ROI = state.imageProc.currentImage(state.imageProc.ROIPositionVector(1):...
        state.imageProc.ROIPositionVector(2),state.imageProc.ROIPositionVector(3):state.imageProc.ROIPositionVector(4),...
        state.imageProc.currentFrame:state.imageProc.numberOfFrames);
    loadImageFromArray('state.imageProc.ROI');
    out=1;
end

if nargout==1
    varargout{1}=out;
end
