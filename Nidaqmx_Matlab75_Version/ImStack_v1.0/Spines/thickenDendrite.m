function thickenDendrite
global gh state

state.imageProc.spine.bwDendrite = bwmorph(state.imageProc.spine.bwSkeleton, 'dilate', state.imageProc.spine.dendriteThickness);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwDendrite);