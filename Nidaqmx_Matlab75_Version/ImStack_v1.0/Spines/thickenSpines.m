function thickenSpines
global gh state

state.imageProc.spine.fatspinesOnly=bwmorph(state.imageProc.spine.spinesOnly, 'thicken', state.imageProc.spine.spineThickness);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.fatspinesOnly);
state.imageProc.spine.fatspinesOnly = double(state.imageProc.spine.fatspinesOnly).*state.imageProc.spine.bw;
state.imageProc.spine.fatspinesOnly = bwmorph(state.imageProc.spine.fatspinesOnly, 'majority',2);
state.imageProc.spine.fatspinesOnly = bwmorph(state.imageProc.spine.fatspinesOnly, 'clean');

set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.fatspinesOnly);

