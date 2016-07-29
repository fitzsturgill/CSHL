function openSpineDataFile(filename)
global gh state

if nargin < 1
	[fname, pname] = uigetfile('*.txt', 'Choose Data to Open...');
	if fname < 1
		return
	else
		filename = [pname fname];
	end	
else
	[pname,fname,ext] = fileparts(filename);
	fname = [fname ext];
end
set([gh.spineDataGUI.Analysis gh.spineDataGUI.add], 'Enable', 'on')
set(gh.spineDataGUI.figure1, 'pointer', 'watch');
cd(pname);
state.imageProc.spineData.spineDataName = fname(1:(end-4));
updateGUIByGlobal('state.imageProc.spineData.spineDataName');

[state.imageProc.spineData.allText, meanDensity, meanDenErr, ...
		overallDensity, state.imageProc.spineData.overallLength, state.imageProc.spineData.overallLenErr,...
		state.imageProc.spineData.meanLength, state.imageProc.spineData.meanLenErr, state.imageProc.spineData.overallVolume, state.imageProc.spineData.overallVolErr, ...
		state.imageProc.spineData.meanVolume, state.imageProc.spineData.meanVolErr, state.imageProc.spineData.HistogramData, state.imageProc.spineData.volumeHist] = readSpineData(filename);


state.imageProc.spineData.overallDensity = overallDensity(1);
state.imageProc.spineData.overallDenErr = 0;
state.imageProc.spineData.meanDenErr = meanDenErr(1);
state.imageProc.spineData.meanDensity = meanDensity(1);
state.imageProc.spineData.overallDensity3 = overallDensity(2);
state.imageProc.spineData.overallDenErr3=  0;
state.imageProc.spineData.meanDenErr3 = meanDenErr(2);
state.imageProc.spineData.meanDensity3 = meanDensity(2);


updateGUIByGlobal('state.imageProc.spineData.meanDensity');
updateGUIByGlobal('state.imageProc.spineData.meanDenErr');
updateGUIByGlobal('state.imageProc.spineData.overallDensity');
updateGUIByGlobal('state.imageProc.spineData.overallDenErr');
updateGUIByGlobal('state.imageProc.spineData.meanDensity3');
updateGUIByGlobal('state.imageProc.spineData.meanDenErr3');
updateGUIByGlobal('state.imageProc.spineData.overallDensity3');
updateGUIByGlobal('state.imageProc.spineData.overallDenErr3');
updateGUIByGlobal('state.imageProc.spineData.meanLength');
updateGUIByGlobal('state.imageProc.spineData.meanLenErr');
updateGUIByGlobal('state.imageProc.spineData.overallLength');
updateGUIByGlobal('state.imageProc.spineData.overallLenErr');
updateGUIByGlobal('state.imageProc.spineData.meanVolume');
updateGUIByGlobal('state.imageProc.spineData.meanVolErr');
updateGUIByGlobal('state.imageProc.spineData.overallVolume');
updateGUIByGlobal('state.imageProc.spineData.overallVolErr');

set(gh.spineDataGUI.figure1, 'pointer', 'arrow');
