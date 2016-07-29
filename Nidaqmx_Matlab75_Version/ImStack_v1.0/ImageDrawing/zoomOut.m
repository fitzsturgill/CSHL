function zoomOut(axis)
global gh state

set(axis, 'XLim', state.imageProc.spine.XLim, 'YLim', state.imageProc.spine.YLim);
