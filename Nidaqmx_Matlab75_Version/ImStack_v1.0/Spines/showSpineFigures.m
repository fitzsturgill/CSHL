function showSpineFigures(name)
global gh state

switch name
	case 'top'
	set(gh.spineGUI.initFigure, 'visible', 'on');
	case 'bot'
	set(gh.spineGUI.initFigure2, 'visible', 'on');
	case 'main'
	set(gh.spineGUI.mainFigure, 'visible', 'on');
case 'preview'
	set(gh.spineGUI.previewFigure, 'visible', 'on');
end
