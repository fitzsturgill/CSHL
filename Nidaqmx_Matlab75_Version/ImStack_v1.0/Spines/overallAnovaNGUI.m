function overallAnovaNGUI(type)
global state gh

p = doOverallAnovaN(type);
if ~strcmp(type, 'all')
	state.imageProc.spineData.annovaPCon = p(1);
	state.imageProc.spineData.annovaP = p(2);
	updateGUIByGlobal('state.imageProc.spineData.annovaPCon');
	updateGUIByGlobal('state.imageProc.spineData.annovaP');
	if length(p) > 2
		state.imageProc.spineData.annovaPIntxn = p(3);
		updateGUIByGlobal('state.imageProc.spineData.annovaPIntxn');
	end
else
	p = mean(p,2);
	state.imageProc.spineData.annovaPCon = p(1);
	state.imageProc.spineData.annovaP = p(2);
	updateGUIByGlobal('state.imageProc.spineData.annovaPCon');
	updateGUIByGlobal('state.imageProc.spineData.annovaP');
	if length(p) > 2
		state.imageProc.spineData.annovaPIntxn = p(3);
		updateGUIByGlobal('state.imageProc.spineData.annovaPIntxn');
	end	
end
