function out = getDataFromExcelSheet(rowStart, rowEnd, columnStart, columnEnd)
global gh state
eval(['out = ddereq(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(rowStart) 'c' num2str(columnStart) ':r' num2str(rowEnd) 'c' num2str(columnEnd)] '''' ');']);
eval(['ddeadv(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(rowStart) 'c' num2str(columnStart) ':r' num2str(rowEnd) 'c' num2str(columnEnd)] '''' ', ''size(state.imageProc.excellink.currentData);'', ''state.imageProc.excellink.currentData'');']);


