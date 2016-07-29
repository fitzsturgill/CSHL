function connectToExcel(filename)
global gh state

try
	state.imageProc.excellink.excelChannel = ddeinit('excel', filename);
catch
	disp('Not a valid File name');
end
