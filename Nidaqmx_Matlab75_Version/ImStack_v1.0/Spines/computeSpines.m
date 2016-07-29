function computeSpines(display)
global gh state

if nargin < 1
	display = 0;
end

[labeled,numObjects] = bwlabel(state.imageProc.spine.fatspinesOnly,8);% Label components.
state.imageProc.spine.numberOfSpines = numObjects;
updateGUIByGlobal('state.imageProc.spine.numberOfSpines');
state.imageProc.spine.spineDensity = state.imageProc.spine.numberOfSpines/state.imageProc.spine.dendriteLength;
updateGUIByGlobal('state.imageProc.spine.spineDensity');

set(state.imageProc.spine.bwimagehandle, 'CData', state.imageProc.spine.fatspinesOnly);


if display == 1
	h=figure('Position', [443   436   570   472]);
	map = [0 0 0; 1 1 1;jet(numObjects+1);1 1 1]; % Create a colormap.
	subplot(1,2,1)
	imshow(state.imageProc.spine.bw+labeled+1,map); 	 % Offset indices to colormap by 1.
	subplot(1,2,2)
	imshow(state.imageProc.spine.fatspinesOnly);
	colormap(map);
	
else
	h=figure;
	map = [0 0 0; 1 1 1;jet(numObjects+1);1 1 1]; % Create a colormap.
	imshow(state.imageProc.spine.bw+labeled+1,map); 
	colormap(map);
	set(gcf,'Position', [443   436   570   472]);
end

state.imageProc.spine.spineFeatures = imfeature(labeled,'all');
