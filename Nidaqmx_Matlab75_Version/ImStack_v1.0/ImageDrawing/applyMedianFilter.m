function applyMedianFilter
global gh state

value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
state.imageProc.filteredImage = [];
sizeArray = size(state.imageProc.currentImage(y:y1,x:x1,state.imageProc.currentFrame:state.imageProc.numberOfFrames),3);
h = waitbar(0,'Filtering of Multi-Tif', 'Name', 'Median Filtering');
total=state.imageProc.numberOfFrames-state.imageProc.currentFrame;
counter=0;
for i = state.imageProc.currentFrame:state.imageProc.numberOfFrames
	waitbar(counter/total,h,['Filtering frame ' num2str(counter) ' of ' num2str(total)]);
	counter=counter+1;
	state.imageProc.filteredImage(:,:,i) = medfilt2(state.imageProc.currentImage(y:y1,x:x1,i),[3 3]);
end
waitbar(1,h,['Loading image']);
loadImageFromArray('state.imageProc.filteredImage');
close(h);