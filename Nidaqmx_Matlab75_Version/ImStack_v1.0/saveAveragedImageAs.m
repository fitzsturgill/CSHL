function saveAveragedImageAs
global gh state

[fname pname] = uiputfile('*.tif', 'Save Averaged Image As...');
% uiputfile has no extension associated with it...

newfilename = [pname fname]; % no extension

value = get(gh.averagingGUI.fileName, 'Value');
image = state.imageProc.cell.currentImage{value};

newImage = genericBin(image, state.imageProc.internal.binX, state.imageProc.internal.binY, state.imageProc.internal.binZ);
% write the file to disk
arrayToTiff(newImage, newfilename, '');
