function loadAveragedImage
global gh state

value = get(gh.averagingGUI.fileName, 'Value');
image = state.imageProc.cell.currentImage{value};

eval(['state.imageProc.internal.averageImage' num2str(state.imageProc.internal.averageImageCounter) ...
	' = genericBin(image, state.imageProc.internal.binX, state.imageProc.internal.binY, state.imageProc.internal.binZ);'])

loadImageFromArray(['state.imageProc.internal.averageImage' num2str(state.imageProc.internal.averageImageCounter)]);
state.imageProc.internal.averageImageCounter = state.imageProc.internal.averageImageCounter+1;