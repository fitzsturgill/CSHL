function reloadSpineImage
global gh state
if state.imageProc.spine.topImage
	
	state.imageProc.spine.maxFlag = 0;
	columns = size(state.imageProc.spine.initialImage,2);
	rows = size(state.imageProc.spine.initialImage,1);
	state.imageProc.spine.totalSpineFrames = size(state.imageProc.spine.initialImage,3);
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames;
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames');
	state.imageProc.spine.currentSpineFrame = state.imageProc.spine.reloadFrameTop;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.currentSpineFrame');
	updateGUIByGlobal('state.imageProc.spine.startMax');
	set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames, 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames 1/state.imageProc.spine.totalSpineFrames]);
	
% 	state.imageProc.spine.lowPixelValue = min(min(min(state.imageProc.spine.initialImage)));
% 	updateGUIByGlobal('state.imageProc.spine.lowPixelValue');
% 	state.imageProc.spine.highPixelValue = .4*double(max(max(max(state.imageProc.spine.initialImage))));
% 	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
	
	set(gh.spineGUI.figure1,'Colormap', makeColormap('gray',8));
	set(gh.spineGUI.initialaxis,'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	state.imageProc.spine.imagehandle = image('CData', state.imageProc.spine.initialImage(:,:,state.imageProc.spine.currentSpineFrame), ...
		'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis, 'ButtonDownFcn', 'spineImageOverFcn');
	state.imageProc.spine.ZCurrentTop = state.imageProc.spine.ZStartTop + (state.imageProc.spine.currentSpineFrame-1)*state.imageProc.spine.zStepSizeTop;
	updateGUIByGlobal('state.imageProc.spine.ZCurrentTop');
	
elseif state.imageProc.spine.bottomImage
	
	state.imageProc.spine.maxFlag2 = 0;
	columns = size(state.imageProc.spine.initialImage2,2);
	rows = size(state.imageProc.spine.initialImage2,1);
	state.imageProc.spine.totalSpineFrames2 = size(state.imageProc.spine.initialImage2,3);
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames2;
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames2');
	state.imageProc.spine.currentSpineFrame2 = state.imageProc.spine.reloadFrameBot;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.currentSpineFrame2');
	updateGUIByGlobal('state.imageProc.spine.startMax');
	set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider2], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames2, 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames2 1/state.imageProc.spine.totalSpineFrames2]);
	
	state.imageProc.spine.lowPixelValue = min(min(min(state.imageProc.spine.initialImage2)));
	updateGUIByGlobal('state.imageProc.spine.lowPixelValue');
	state.imageProc.spine.highPixelValue = .1*double(max(max(max(state.imageProc.spine.initialImage2))));
	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
	
	set(gh.spineGUI.figure1,'Colormap', makeColormap('gray',8));
	set(gh.spineGUI.initialaxis2,'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	state.imageProc.spine.imagehandle2 = image('CData', state.imageProc.spine.initialImage2(:,:,state.imageProc.spine.currentSpineFrame2), ...
		'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis2, 'ButtonDownFcn', 'spineImageOverFcn');
	state.imageProc.spine.ZCurrentBot = state.imageProc.spine.ZStartBot + (state.imageProc.spine.currentSpineFrame2-1)*state.imageProc.spine.zStepSizeBot;
	updateGUIByGlobal('state.imageProc.spine.ZCurrentBot');
else
	return
end

