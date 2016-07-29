function showInitialSpines
global gh state

state.imageProc.spine.bwDendrite=~state.imageProc.spine.bwDendrite;
state.imageProc.spine.bwDendrite=double(state.imageProc.spine.bwDendrite);

state.imageProc.spine.bw=double(state.imageProc.spine.bw);
state.imageProc.spine.spinesOnly = state.imageProc.spine.bw.*state.imageProc.spine.bwDendrite;
state.imageProc.spine.spinesOnly = bwmorph(state.imageProc.spine.spinesOnly, 'majority',2);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.spinesOnly);
