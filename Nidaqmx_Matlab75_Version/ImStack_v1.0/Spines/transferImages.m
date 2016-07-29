function transferImages
	
	global gh state
	
%	set(gh.spineGUI.mainFigure, 'visible', 'on');
	if state.imageProc.spine.auto
		showBWImage;
		return
	end
	if state.imageProc.spine.eraseSpines
		try
			removeLines;
		end
	end
	
	if state.imageProc.spine.topImage
	%	state.imageProc.spine.baseName=[state.imageProc.spine.loadedFileNameTop '0'];
		updateGUIByGlobal('state.imageProc.spine.baseName');
		if state.imageProc.spine.maxFlag == 0 % max not dsipalyed
			[x,x1,y,y1] = getCurrentAxisLimits(gh.spineGUI.initialaxis);
			set(gh.spineGUI.mainAxes, 'YLim', [y y1], 'XLim', [x x1], 'Ydir', 'reverse', 'Clim', [state.imageProc.spine.lowPixelValue ...
					state.imageProc.spine.highPixelValue]);
			state.imageProc.spine.maniImgData = state.imageProc.spine.initialImage(y:y1,x:x1,state.imageProc.spine.currentSpineFrame);
			set(state.imageProc.spine.mainImage, 'CData', state.imageProc.spine.maniImgData, 'YData', [y y1], 'XData', [x x1]);
			state.imageProc.spine.from = state.imageProc.spine.ZCurrentTop;
			state.imageProc.spine.to = state.imageProc.spine.ZCurrentTop;
		else
			[x,x1,y,y1] = getCurrentAxisLimits(gh.spineGUI.initialaxis);
			set(gh.spineGUI.mainAxes, 'YLim', [y y1], 'XLim', [x x1], 'Ydir', 'reverse', 'Clim', [state.imageProc.spine.lowPixelValue ...
					state.imageProc.spine.highPixelValue]);
			state.imageProc.spine.maniImgData = state.imageProc.spine.maxProjection(y:y1,x:x1);
			set(state.imageProc.spine.mainImage, 'CData', state.imageProc.spine.maniImgData, 'YData', [y y1], 'XData', [x x1]);
			state.imageProc.spine.from = state.imageProc.spine.from1;
			state.imageProc.spine.to = state.imageProc.spine.to1;
		end
		
	elseif state.imageProc.spine.bottomImage
	%	state.imageProc.spine.baseName=[state.imageProc.spine.loadedFileNameBot '0'];
		updateGUIByGlobal('state.imageProc.spine.baseName');
		
		if state.imageProc.spine.maxFlag2 == 0 % max not dsipalyed
			[x,x1,y,y1] = getCurrentAxisLimits(gh.spineGUI.initialaxis2);
			set(gh.spineGUI.mainAxes, 'YLim', [y y1], 'XLim', [x x1], 'Ydir', 'reverse', 'Clim', [state.imageProc.spine.lowPixelValue ...
					state.imageProc.spine.highPixelValue]);
			state.imageProc.spine.maniImgData =  state.imageProc.spine.initialImage2(y:y1,x:x1,state.imageProc.spine.currentSpineFrame2);
			state.imageProc.spine.mainImage = image('CData',state.imageProc.spine.maniImgData,...
				'CDataMapping', 'scaled','Parent', gh.spineGUI.mainAxes,'YData', [y y1], 'XData', [x x1]);
			state.imageProc.spine.from = state.imageProc.spine.ZCurrentBot;
			state.imageProc.spine.to = state.imageProc.spine.ZCurrentBot;
		else
			[x,x1,y,y1] = getCurrentAxisLimits(gh.spineGUI.initialaxis2);
			set(gh.spineGUI.mainAxes, 'YLim', [y y1], 'XLim', [x x1], 'Ydir', 'reverse', 'Clim', [state.imageProc.spine.lowPixelValue ...
					state.imageProc.spine.highPixelValue]);
			state.imageProc.spine.maniImgData = state.imageProc.spine.maxProjection2(y:y1,x:x1);
			state.imageProc.spine.mainImage = image('CData',state.imageProc.spine.maniImgData ,...
				'CDataMapping', 'scaled','Parent', gh.spineGUI.mainAxes,'YData', [y y1], 'XData', [x x1]);
			state.imageProc.spine.from = state.imageProc.spine.from2;
			state.imageProc.spine.to = state.imageProc.spine.to2;
		end
	else
		display('No Image Loaded');
		return
	end