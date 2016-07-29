function applyWienerFilter
global gh state

value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
state.imageProc.filteredImage = [];
sizeArray = size(state.imageProc.currentImage(y:y1,x:x1,state.imageProc.currentFrame:state.imageProc.numberOfFrames),3);
for i = state.imageProc.currentFrame:state.imageProc.numberOfFrames
	state.imageProc.filteredImage(:,:,i) = wiener2(state.imageProc.currentImage(y:y1,x:x1,i),[5 5]);
end

loadImageFromArray('state.imageProc.filteredImage');
