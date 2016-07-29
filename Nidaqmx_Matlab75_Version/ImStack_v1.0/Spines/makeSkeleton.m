function makeSkeleton(image)
global state gh

h = waitbar(0,'Calculating Medial Projection..', 'Name', 'Spine Analysis Status');

waitbar(1/8,h, 'Cleaning Stray Pixels');
state.imageProc.spine.bwSkeleton = bwmorph(image, 'clean');
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;
waitbar(2/8,h, 'Majority Editing');
state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'majority');
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;
waitbar(3/8,h, 'Initial Skeleton');
state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'skel', inf);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;
state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'fill', 5);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;
state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'thin', 5);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;
waitbar(4/8,h, 'Removing Spurs');
state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'spur', inf);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;
waitbar(5/8,h, 'Cleaning Stray Pixels Again');
state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'clean');
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;
waitbar(6/8,h, 'Connecting Regions');
state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'fill', 5);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;
waitbar(7/8,h, 'Thinning Dendrite');
state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'thin', 5);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;
%state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'shrink',1);
state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'clean',5);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);
drawnow;

state.imageProc.spine.bwSkeleton = (double(state.imageProc.spine.bwSkeleton));
waitbar(8/8,h, 'Updating Dendrite Statistics');

state.imageProc.spine.dendriteLength = bwarea(state.imageProc.spine.bwSkeleton)*(state.imageProc.spine.micronsperpixelX+state.imageProc.spine.micronsperpixelY)/2;
updateGUIByGlobal('state.imageProc.spine.dendriteLength');

state.imageProc.spine.bwSkeleton = bwmorph(state.imageProc.spine.bwSkeleton, 'fill', 5);
set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.bwSkeleton);

close(h);