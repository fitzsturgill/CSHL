function loadSpineImageBot
global gh state

% Call load image ui
try
	cd(state.imageProc.spine.openPath);
end

[fname, pname] = uigetfile({'*.tif;*.CFD'}, 'Choose image to load to bottom axis');

	if ~isnumeric(fname)
		filename = [pname fname];
	else
		return
	end

	[path,name,ext] = fileparts(filename);
	
	
	switch ext 
		case '.tif'
			try
				state.imageProc.spine.headerBot = readImageHeaderTif(filename);
				state.imageProc.spine.zStepSizeBot = valueFromHeaderString('state.acq.zStepSize', state.imageProc.spine.headerBot);
				state.imageProc.spine.ZStartBot = valueFromHeaderString('state.motor.relZPosition', state.imageProc.spine.headerBot);
				state.imageProc.spine.ZCurrentBot = state.imageProc.spine.ZStartBot;
				updateGUIByGlobal('state.imageProc.spine.ZCurrentBot');
			catch
				state.imageProc.spine.headerBot = '';
				state.imageProc.spine.zStepSizeBot = -1;
				state.imageProc.spine.ZStartBot = 0;
				state.imageProc.spine.ZCurrentBot = state.imageProc.spine.ZStartBot;
				updateGUIByGlobal('state.imageProc.spine.ZCurrentBot');
			end
			state.imageProc.spine.initialImage2 = opentif(filename);
		case '.CFD'
			state.imageProc.spine.headerBot = '';
			state.imageProc.spine.zStepSizeBot = -1;
			state.imageProc.spine.ZStartBot = 0;
			state.imageProc.spine.ZCurrentBot = state.imageProc.spine.ZStartBot;
			updateGUIByGlobal('state.imageProc.spine.ZCurrentBot');
  		    state.imageProc.spine.initialImage2 = openACFDSpine([pname fname],state.imageProc.cfd.numberofChannels);
		case '.cfd'
			state.imageProc.spine.headerBot = '';
			state.imageProc.spine.zStepSizeBot = -1;
			state.imageProc.spine.ZStartBot = 0;
			state.imageProc.spine.ZCurrentBot = state.imageProc.spine.ZStartBot;
			updateGUIByGlobal('state.imageProc.spine.ZCurrentBot');
		   state.imageProc.spine.initialImage2 = openACFDSpine([pname fname],state.imageProc.cfd.numberofChannels);
		otherwise
			disp('File Not Recognized');
			return
	end
	state.imageProc.spine.openPath = pname;
	state.imageProc.spine.loadedFileNameBot = name;
	state.imageProc.spine.maxFlag2 = 0;
	state.imageProc.spine.topImage = 0;
	state.imageProc.spine.bottomImage = 1;
	updateGUIByGlobal('state.imageProc.spine.topImage');
	updateGUIByGlobal('state.imageProc.spine.bottomImage');
	set(gh.spineGUI.bottomImage, 'Enable', 'on');
	
	columns = size(state.imageProc.spine.initialImage2,2);
	rows = size(state.imageProc.spine.initialImage2,1);
	state.imageProc.spine.totalSpineFrames2 = size(state.imageProc.spine.initialImage2,3);
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames2;
	state.imageProc.spine.currentSpineFrame2 = 1;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.currentSpineFrame2');
	updateGUIByGlobal('state.imageProc.spine.startMax');
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames2');
	if  state.imageProc.spine.totalSpineFrames2 > 1
		set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider2], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames2, 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames2 1/state.imageProc.spine.totalSpineFrames2]);
	else
		set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider2], 'Min', 1, 'Max', 1.001, 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames2 1/state.imageProc.spine.totalSpineFrames2]);
	end
	
	
	state.imageProc.spine.lowPixelValue = min(min(min(state.imageProc.spine.initialImage2)));
	updateGUIByGlobal('state.imageProc.spine.lowPixelValue');
	state.imageProc.spine.highPixelValue = .1*double(max(max(max(state.imageProc.spine.initialImage2))));
	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
	
	set(gh.spineGUI.initialaxis2,'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue], 'Parent', gh.spineGUI.initFigure2);
	set(gh.spineGUI.initFigure2, 'visible','on');
	state.imageProc.spine.XLim2 = [1 columns];
	state.imageProc.spine.YLim2 = [1 rows];
	state.imageProc.spine.spineThreshold = state.imageProc.spine.highPixelValue;
	updateGUIByGlobal('state.imageProc.spine.spineThreshold');
	state.imageProc.spine.imagehandle2 = image('CData', state.imageProc.spine.initialImage2(:,:,1), ...
		'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis2, 'ButtonDownFcn', 'spineImageOverFcnBot');
	
	
	
