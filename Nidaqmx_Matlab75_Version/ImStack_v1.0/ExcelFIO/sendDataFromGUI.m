function sendDataFromGUI
global gh state

eval(['state.imageProc.excellink.currentData = ' state.imageProc.excellink.dataToWrite ';']);
eval(['sendDataToExcelSheet(' state.imageProc.excellink.dataToWrite ');']);