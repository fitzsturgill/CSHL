function addCurrentImgToStack
global state
% THis function will read the current frame into a stack that can be loaded as an image in the image prcessing
% environment.

state.imageProc.imgStack = cat(3,state.imageProc.imgStack,state.imageProc.currentImage(:,:,...
	state.imageProc.currentFrame));
