function open3DMAReader
global state gh

seeGUI('gh.spineDataGUI.figure1');
set([gh.spineDataGUI.Analysis gh.spineDataGUI.add], 'Enable', 'off');
set(gh.spineDataGUI.combSpData, 'Enable', 'off');
hideGUI('gh.imageParsingGUI.figure1');
