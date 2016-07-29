function openSpineAnalysis
global gh state

seeGUI('gh.spineGUI.figure1');
set([gh.spineGUI.initFigure gh.spineGUI.mainFigure], 'Visible', 'on');
state.imageProc.cfd.numberofChannels =1;
updateGUIByGlobal('state.imageProc.cfd.numberofChannels');

button = questdlg('Do you want to open an Excel Sheet to store spine data?',...
'Open Excel File?','Yes','No','Yes');

if strcmp(button,'Yes')
   loadExcelFile;
elseif strcmp(button,'No')
   return
end

