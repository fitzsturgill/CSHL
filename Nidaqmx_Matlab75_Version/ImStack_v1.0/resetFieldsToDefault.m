function resetFieldsToDefault
global gh state

	state.imageProc.currentImage= [];
	state.imageProc.fileName='Load Image...';
	state.imageProc.numberOfFrames=1;
	updateGUIByGlobal('state.imageProc.numberOfFrames');
	state.imageProc.parsing.pixelsPerLine=256;
	updateGUIByGlobal('state.imageProc.parsing.pixelsPerLine');
	state.imageProc.parsing.linesPerFrame=256;
	updateGUIByGlobal('state.imageProc.parsing.linesPerFrame');
	state.imageProc.montageStart=1;
	updateGUIByGlobal('state.imageProc.montageStart');
	state.imageProc.montageEnd=1;
	updateGUIByGlobal('state.imageProc.montageEnd');
	state.imageProc.parsing.numberOfChannels=1;
	updateGUIByGlobal('state.imageProc.parsing.numberOfChannels');
	state.imageProc.parsing.numberOfFrames=1;
	updateGUIByGlobal('state.imageProc.parsing.numberOfFrames');
	state.imageProc.parsing.scanRotation=0;
	updateGUIByGlobal('state.imageProc.parsing.scanRotation');
	state.imageProc.parsing.averaged=0;
	updateGUIByGlobal('state.imageProc.parsing.averaged');
	state.imageProc.parsing.numberOfZSlices=1;
	updateGUIByGlobal('state.imageProc.parsing.numberOfZSlices');
	state.imageProc.maxEnd=1;
	updateGUIByGlobal('state.imageProc.maxEnd');
	state.imageProc.maxStart=1;
	updateGUIByGlobal('state.imageProc.maxStart');
	state.imageProc.movieEnd=1;
	updateGUIByGlobal('state.imageProc.movieEnd');
	state.imageProc.movieStart=1;
	updateGUIByGlobal('state.imageProc.movieStart');
	state.imageProc.currentFrame=1;
	updateGUIByGlobal('state.imageProc.currentFrame');
	state.imageProc.parsing.header='Image Description';
	updateGUIByGlobal('state.imageProc.parsing.header');
	value = get(gh.montageGUI.fileName, 'Value');
	
	
	set(gh.montageGUI.montageStartSlider, 'Min', 1 , 'Max',  1001, 'SliderStep', [.001 .001]);
	set(gh.montageGUI.montageEndSlider, 'Min', 1 , 'Max',  1001, 'SliderStep', [.001 .001]);

	set(gh.maxProjectionGUI.maxStartSlider, 'Min', 1 , 'Max',  1001, 'SliderStep', [.001 .001]);
	set(gh.maxProjectionGUI.maxEndSlider, 'Min', 1 , 'Max',  1001, 'SliderStep', [.001 .001]);

	set(gh.movieGUI.movieStartSlider , 'Min', 1 , 'Max',  1001, 'SliderStep', [.001 .001]);
	set(gh.movieGUI.movieEndSlider , 'Min', 1 , 'Max',  1001, 'SliderStep', [.001 .001]);

	set(gh.imageProcessingGUI.currentFrameSlider , 'Min', 1 , 'Max',  1001, 'SliderStep', [.001 .001]);
	set(gh.imageProcessingGUI.totalFramesSlider, 'Min', 1 , 'Max',  1001, 'SliderStep', [.001 .001]);
	
	state.imageProc.internal.imageCounter = 1;
		
		