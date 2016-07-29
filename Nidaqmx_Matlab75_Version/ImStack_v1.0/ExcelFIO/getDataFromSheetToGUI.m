function getDataFromSheetToGUI
global gh state

state.imageProc.excellink.currentData = getDataFromExcelSheet(state.imageProc.excellink.startRowRead, ...
	state.imageProc.excellink.endRowRead, state.imageProc.excellink.startColRead, state.imageProc.excellink.endColRead);

evalin('base', [state.imageProc.excellink.readDataArray '=state.imageProc.excellink.currentData;']);
disp(['Made variable ' state.imageProc.excellink.readDataArray ' from excel sheet']);

