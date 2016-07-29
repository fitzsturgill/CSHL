function sendDataToExcelSheet(data)
global gh state

rows = size(data,1);
columns = size(data,2);

eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r1c1:r' num2str(rows) 'c' num2str(columns)] '''' ', data);']);
eval(['ddeadv(state.imageProc.excellink.excelChannel, ' '''' ['r1c1:r' num2str(rows) 'c' num2str(columns) ] '''' ', ''size(state.imageProc.excellink.currentData)'', ''state.imageProc.excellink.currentData'');']);

	