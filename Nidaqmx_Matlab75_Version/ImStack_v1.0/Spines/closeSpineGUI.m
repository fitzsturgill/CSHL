function closeSpineGUI
global gh state


set([gh.spineGUI.initFigure gh.spineGUI.initFigure2 gh.spineGUI.mainFigure gh.spineGUI.previewFigure], 'visible', 'off');
hideGUI('gh.spineGUI.figure1');