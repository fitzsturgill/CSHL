function openAndLoadSpineData
global gh state

if isempty(state.imageProc.spineData.openPath)
	curDirSpineData;
else
	cd(state.imageProc.spineData.openPath);
end

filenames = selectFilesFromList(state.imageProc.spineData.openPath, '.txt');

if isempty(filenames)
	return
end

h = waitbar(0,'Starting To Load Spine Data Files....', 'Name', '3DMA Spine Data File Reader', ...
	'Position', [164.250000000000   693.750000000000   270.000000000000   056.250000000000], ...
	'Pointer', 'watch');
total = length(filenames);

for i = 1:total
	waitbar((i-1)/total,h, filenames{i});
	openSpineDataFile([state.imageProc.spineData.openPath filenames{i}]);
	waitbar(i/total,h);
	addSpineDataToList;
end

close(h);
	
