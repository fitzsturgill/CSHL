function loadSpineImage
global gh state

% Call load image ui
% Call load image ui
%	saveMainAxis;
disp('********** autosave picture off **********');

	try
		cd(state.imageProc.spine.openPath);
	end
	[fname, pname] = uigetfile({[state.imageProc.spine.baseName '*max*.tif;']}, 'Choose image to load to top axis');
%	[fname, pname] = uigetfile({[state.imageProc.spine.baseName '*.tif;']}, 'Choose image to load to top axis');

	if ~isnumeric(fname)
		filename = [pname fname];
	else
		return
	end

	[path,name,ext] = fileparts(filename);
	
	switch ext 
		case '.tif'
			try
				state.imageProc.spine.headerTop = readImageHeaderTif(filename);
				state.imageProc.spine.zStepSizeTop = valueFromHeaderString('state.acq.zStepSize', state.imageProc.spine.headerTop);
				state.imageProc.spine.ZStartTop = valueFromHeaderString('state.motor.relZPosition', state.imageProc.spine.headerTop);
				state.imageProc.spine.ZCurrentTop = state.imageProc.spine.ZStartTop;
				updateGUIByGlobal('state.imageProc.spine.ZCurrentTop');
			catch
				state.imageProc.spine.headerTop = '';
				state.imageProc.spine.zStepSizeTop = -1;
				state.imageProc.spine.ZStartTop = 0;
				state.imageProc.spine.ZCurrentTop = state.imageProc.spine.ZStartTop;
				updateGUIByGlobal('state.imageProc.spine.ZCurrentTop');
			end
            chanOn=[0 0 0];
            chanOn(1)=valueFromHeaderString('state.acq.savingChannel1', state.imageProc.spine.headerTop);
            chanOn(2)=valueFromHeaderString('state.acq.savingChannel2', state.imageProc.spine.headerTop);
            chanOn(3)=valueFromHeaderString('state.acq.savingChannel3', state.imageProc.spine.headerTop);
            nChan=length(find(chanOn));
            if nChan>1
    			initialImage = opentif(filename);
 	            nSlices=size(initialImage,3)/nChan;
                frames=nChan*[1:nSlices]-1;
                state.imageProc.spine.initialImage=initialImage(:,:,frames);
            else
    			state.imageProc.spine.initialImage = opentif(filename);
            end
			state.imageProc.spine.from=0;
			state.imageProc.spine.to=0;
			
		case '.CFD'
			state.imageProc.spine.headerTop = '';
			state.imageProc.spine.zStepSizeTop = -1;
			state.imageProc.spine.ZStartTop = 0;
			state.imageProc.spine.ZCurrentTop = state.imageProc.spine.ZStartTop;
			updateGUIByGlobal('state.imageProc.spine.ZCurrentTop');
			[pname fname]
  		    state.imageProc.spine.initialImage = openACFDSpine([pname fname],state.imageProc.cfd.numberofChannels);
		case '.cfd'
			state.imageProc.spine.headerTop = '';
			state.imageProc.spine.zStepSizeTop = -1;
			state.imageProc.spine.ZStartTop = 0;
			state.imageProc.spine.ZCurrentTop = state.imageProc.spine.ZStartTop;
			updateGUIByGlobal('state.imageProc.spine.ZCurrentTop');
  		    state.imageProc.spine.initialImage = openACFDSpine([pname fname],state.imageProc.cfd.numberofChannels);
		otherwise
			disp('File Not Recognized');
			return
	end
	
	state.imageProc.spine.openPath = pname;
	state.imageProc.spine.loadedFileNameTop = name;
	state.imageProc.spine.topImage = 1;
	state.imageProc.spine.bottomImage = 0;
	updateGUIByGlobal('state.imageProc.spine.topImage');
	updateGUIByGlobal('state.imageProc.spine.bottomImage');
	set(gh.spineGUI.topImage, 'Enable', 'on');
	state.imageProc.spine.maxFlag = 0;
	columns = size(state.imageProc.spine.initialImage,2);
	rows = size(state.imageProc.spine.initialImage,1);
	state.imageProc.spine.totalSpineFrames = size(state.imageProc.spine.initialImage,3);
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames;
	state.imageProc.spine.currentSpineFrame = 1;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.currentSpineFrame');
	updateGUIByGlobal('state.imageProc.spine.startMax');
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames');
	if  state.imageProc.spine.totalSpineFrames > 1
		set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames, 'SliderStep', ...
			[1/state.imageProc.spine.totalSpineFrames 1/state.imageProc.spine.totalSpineFrames]);
	else
		set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider], 'Min', 1, 'Max', 1.001, 'SliderStep', ...
			[1/state.imageProc.spine.totalSpineFrames 1/state.imageProc.spine.totalSpineFrames]);
	end
	
% 	state.imageProc.spine.lowPixelValue = min(min(min(state.imageProc.spine.initialImage)));
% 	updateGUIByGlobal('state.imageProc.spine.lowPixelValue');
% 	state.imageProc.spine.highPixelValue = .1*double(max(max(max(state.imageProc.spine.initialImage))));
% 	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
	
	set(gh.spineGUI.initialaxis,'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue], 'Parent', gh.spineGUI.initFigure);
	set(gh.spineGUI.initFigure, 'visible','on');
	state.imageProc.spine.spineThreshold = state.imageProc.spine.highPixelValue;
	updateGUIByGlobal('state.imageProc.spine.spineThreshold');
	state.imageProc.spine.XLim = [1 columns];
	state.imageProc.spine.YLim = [1 rows];
	set(state.imageProc.spine.imagehandle, 'CData', state.imageProc.spine.initialImage(:,:,1));
% 	state.imageProc.spine.imagehandle = image('CData', state.imageProc.spine.initialImage(:,:,1), ...
% 		'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis, 'ButtonDownFcn', 'spineImageOverFcn');
	
	
	resetSpines;
	resetDendrite;
	incrementExcelRow;
