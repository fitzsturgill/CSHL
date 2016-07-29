function openExcelSheet(filename)
global gh state
% 
% state.imageProc.excellink.Handle = actxserver('Excel.Application');
% set(state.imageProc.excellink.Handle, 'Visible', 1);
% state.imageProc.excellink.Workbook = state.imageProc.excellink.Handle.Workbooks;
% invoke(state.imageProc.excellink.Workbook, 'Open', filename);
% state.imageProc.excellink.Sheets = state.imageProc.excellink.Handle.ActiveWorkBook.Sheets;
% sheet1 = get(state.imageProc.excellink.Sheets , 'Item', 1);
% invoke(sheet1, 'Activate');
% state.imageProc.excellink.Activesheet =  state.imageProc.excellink.Handle.Activesheet;

eval([' ! ' filename '&']); 


