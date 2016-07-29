function updateBaseName
global gh state

state.imageProc.baseName = get(gh.fileCounterGUI.baseName, 'String');
updateGUIByGlobal('state.imageProc.baseName');