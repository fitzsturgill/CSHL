function saveCurrentImageAsTif
global state gh data
% Will save current ROI to file...
try
	cd(state.imageProc.internal.savePath);
catch
end

[fname, pname] = uiputfile({'*.tif'},'Save Region of Interest as...' );
if isnumeric(fname)
	return
else
	filename=[pname fname];
end

state.imageProc.internal.savePath = pname;
updateGUIByGlobal('state.imageProc.internal.savePath');
state.imageProc.internal.saveBaseName = fname;
updateGUIByGlobal('state.imageProc.internal.saveBaseName');
state.imageProc.internal.updateSavePath == 1;

value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
data=state.imageProc.cell.currentImage{value}(y:y1, x:x1,state.imageProc.cell.currentFrame{value}:state.imageProc.cell.numberOfFrames{value});
data=uint16(data);
arrayToTiff(data, filename, '');



