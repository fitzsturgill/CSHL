function annovaGUI(type)
global gh state

p = doMultiAnova(type);
if isempty(p)
	return
end

switch type
case 'density2d'
	state.imageProc.spineData.annova2d = p;
	updateGUIByGlobal('state.imageProc.spineData.annova2d');
case 'density3d'
	state.imageProc.spineData.annova3d = p;
	updateGUIByGlobal('state.imageProc.spineData.annova3d');
case 'length'
	state.imageProc.spineData.annovaLen = p;
	updateGUIByGlobal('state.imageProc.spineData.annovaLen');
case 'volume'
	state.imageProc.spineData.annovaVol = p;
	updateGUIByGlobal('state.imageProc.spineData.annovaVol');
case 'all'
	state.imageProc.spineData.annova2d = p(1);
	state.imageProc.spineData.annova3d = p(2);
	state.imageProc.spineData.annovaLen = p(3);
	state.imageProc.spineData.annovaVol = p(4);
	updateGUIByGlobal('state.imageProc.spineData.annova2d');
	updateGUIByGlobal('state.imageProc.spineData.annova3d');
	updateGUIByGlobal('state.imageProc.spineData.annovaLen');
	updateGUIByGlobal('state.imageProc.spineData.annovaVol');	
end
