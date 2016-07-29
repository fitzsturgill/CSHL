function getArea(axis)
global gh state 

% try
% 	for i = 1:length(state.imageProc.spine.dendriteLines)
% 		set(state.imageProc.spine.dendriteLines(i), 'Visible', 'off');
% 	end
% end

[X,Y] = getptsSpine(axis)
X=[X' X(1)];
Y=[Y' Y(1)];

A=polyarea(X,Y);
state.imageProc.spine.numberOfSpines = state.imageProc.spine.numberOfSpines + 1;
updateGUIByGlobal('state.imageProc.spine.numberOfSpines');
disp(['measured area = ' num2str(A)]);

state.imageProc.spine.spineLengths(end+1) = A;

saveSpineDataToExcel;