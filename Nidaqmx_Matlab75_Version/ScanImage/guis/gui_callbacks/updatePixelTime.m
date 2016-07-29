function updatePixelTime(handle)
global state

state.acq.pixelTime = ((state.acq.fillFraction*state.acq.msPerLine)/state.acq.pixelsPerLine);
updateGUIByGlobal('state.acq.pixelTime');



