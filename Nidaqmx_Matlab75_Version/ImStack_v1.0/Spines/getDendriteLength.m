function getDendriteLength(axis)
global gh state 

% try
% 	for i = 1:length(state.imageProc.spine.dendriteLines)
% 		set(state.imageProc.spine.dendriteLines(i), 'Visible', 'off');
% 	end
% end

XLim = get(axis, 'XLim');
YLim = get(axis, 'YLim');

[X,Y] = getptsSpine(axis);

points = [X Y];
totalLines= length(X)-1;
total = 0;
len=size(state.imageProc.spine.dendriteLines,2);
for i = 1:totalLines
	x = points(i:(i+1),1);
	y = points(i:(i+1),2);
	a = x(2)-x(1);
	b = y(2) - y(1);
	length = ((state.imageProc.spine.micronsperpixelX*a)^2 + (state.imageProc.spine.micronsperpixelY*b)^2)^.5;
	total = total + length;
	if state.imageProc.spine.showDendrite == 0
		state.imageProc.spine.dendriteLines(i+len)=line(x, y, 'tag', 'Dendrite Length' , 'color', [0 1 0],'linewidth', 4, 'visible', 'off');
	else
		state.imageProc.spine.dendriteLines(i+len)=line(x, y, 'tag', 'Dendrite Length' , 'color', [0 1 0],'linewidth', 4, 'visible', 'on');
	end
end
state.imageProc.spine.dendriteLength = state.imageProc.spine.dendriteLength + total;
updateGUIByGlobal('state.imageProc.spine.dendriteLength');

state.imageProc.spine.spineDensity = state.imageProc.spine.numberOfSpines/state.imageProc.spine.dendriteLength;
updateGUIByGlobal('state.imageProc.spine.spineDensity');
drawnow;
pause(0.01);
saveSpineDataToExcel;

	