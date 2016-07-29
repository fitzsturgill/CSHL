function removeLines
global state

for i = 1:length(state.imageProc.spine.linehandles)
	set(state.imageProc.spine.linehandles(i), 'Visible', 'off');
	set(state.imageProc.spine.texthandles(i), 'Visible', 'off');
end
for i = 1:length(state.imageProc.spine.dendriteLines)
	set(state.imageProc.spine.dendriteLines(i), 'Visible', 'off');
end
