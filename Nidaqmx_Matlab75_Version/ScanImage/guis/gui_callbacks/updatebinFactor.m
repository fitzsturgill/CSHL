function updatebinFactor(handle)
global state

state.acq.binFactor = (state.acq.samplesAcquiredPerLine/state.acq.pixelsPerLine);
updateGUIByGlobal('state.acq.binFactor');
