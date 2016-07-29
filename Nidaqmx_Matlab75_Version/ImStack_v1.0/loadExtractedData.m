function loadExtractedData
global state gh

value = get(gh.imageProcessingGUI.fileName, 'Value');

image = state.imageProc.cell.currentImage{value};
numberofchannels = state.imageProc.parsing.cell.numberOfChannels{value};

state.imageProc.internal.newimage = extractData(image, numberofchannels);

for channelcounter = 1:numberofchannels
	loadImageFromArray(['state.imageProc.internal.newimage{' num2str(channelcounter) '}']);
end