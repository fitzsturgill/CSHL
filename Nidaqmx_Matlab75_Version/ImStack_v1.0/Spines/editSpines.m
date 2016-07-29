function editSpines
global gh state

computeSpines(1);
state.imageProc.spine.removedSpines = bwselect(8);
state.imageProc.spine.fatspinesOnly = double(state.imageProc.spine.fatspinesOnly) - double(state.imageProc.spine.removedSpines);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.fatspinesOnly);
close(gcf);
computeSpines;

