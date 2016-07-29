function loadExcelFile
global gh state

[fname, pname] = uigetfile('*.xls', 'Choose Excel File to Open...');
if pname == 0
	return
else
	state.imageProc.excellink.excelFileName = [pname fname];
	updateGUIByGlobal('state.imageProc.excellink.excelFileName');
	openExcelSheet(state.imageProc.excellink.excelFileName);
end

button = questdlg(['Do you want to connect to ' [pname fname] '?'] ,...
'Wait For File To Load Before Connecting!','Yes', 'No','Yes');
if strcmp(button,'Yes')
	connectToExcel(state.imageProc.excellink.excelFileName);
elseif strcmp(button,'No')
   return
end

