function showCurrentFrame(Value)
global state gh
state.imageProc.currentFrame = round(state.imageProc.currentFrame);
set(state.imageProc.internal.imagehandle{Value}, 'CData', state.imageProc.currentImage(:,:,state.imageProc.currentFrame));