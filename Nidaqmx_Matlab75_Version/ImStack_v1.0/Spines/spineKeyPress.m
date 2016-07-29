function spineKeyPress
global gh state
val = double(get(gcbo,'CurrentCharacter'));

try
	switch val
	case 28 %left arrow
		px=get(gh.spineGUI.initialaxis, 'Xlim');
		if px(1)<=1
			beep
			return
		end
		px2=[max(1, px(1)+round((px(2)-px(1))/5)-(px(2)-px(1))) ...
				px(1)+round((px(2)-px(1))/5)];
		py2=get(gh.spineGUI.initialaxis, 'Ylim');
		set(gh.spineGUI.initialaxis, 'Xlim', px2, 'Ylim', py2);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 29	% right arrow
		px=get(gh.spineGUI.initialaxis, 'Xlim');
		if px(2)>=size(state.imageProc.spine.initialImage,2)
			beep
			return
		end
		px2=[px(2)-round((px(2)-px(1))/5) ...
				min(px(2)-round((px(2)-px(1))/5)+(px(2)-px(1)), size(state.imageProc.spine.initialImage,2))]
		py2=get(gh.spineGUI.initialaxis, 'Ylim');
		set(gh.spineGUI.initialaxis, 'Xlim', px2, 'Ylim', py2);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 30	% up arrow
		px2=get(gh.spineGUI.initialaxis, 'Xlim');
		py=get(gh.spineGUI.initialaxis, 'Ylim');
		if py(1)<=1
			beep
			return
		end
		py2=[max(1, py(1)+round((py(2)-py(1))/5)-(py(2)-py(1))) ...
				py(1)+round((py(2)-py(1))/5)];
		set(gh.spineGUI.initialaxis, 'Xlim', px2, 'Ylim', py2);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 31	% down arrow
		px2=get(gh.spineGUI.initialaxis, 'Xlim');
		py=get(gh.spineGUI.initialaxis, 'Ylim');
		if py(2)>=size(state.imageProc.spine.initialImage,1)
			beep
			return
		end
		py2=[py(2)-round((py(2)-py(1))/5) ...
				min(py(2)-round((py(2)-py(1))/5)+(py(2)-py(1)), size(state.imageProc.spine.initialImage,1))];
		set(gh.spineGUI.initialaxis, 'Xlim', px2, 'Ylim', py2);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 33
		px=get(gh.spineGUI.initialaxis, 'Xlim');
		px2=[px(1) px(1)+round(1.1*(px(2)-px(1))/2)];
		py=get(gh.spineGUI.initialaxis, 'Ylim');
		py2=[py(1) py(1)+round(1.1*(py(2)-py(1))/2)];
		set(gh.spineGUI.initialaxis, 'Xlim', px2, 'Ylim', py2);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 64
		px=get(gh.spineGUI.initialaxis, 'Xlim');
		px2=[px(1)+round(0.9*(px(2)-px(1))/2) px(2)];
		py=get(gh.spineGUI.initialaxis, 'Ylim');
		py2=[py(1) py(1)+round(1.1*(py(2)-py(1))/2)];
		set(gh.spineGUI.initialaxis, 'Xlim', px2, 'Ylim', py2);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 35
		px=get(gh.spineGUI.initialaxis, 'Xlim');
		px2=[px(1) px(1)+round(1.1*(px(2)-px(1))/2)];
		py=get(gh.spineGUI.initialaxis, 'Ylim');
		py2=[py(1)+round(0.9*(py(2)-py(1))/2) py(2)];
		set(gh.spineGUI.initialaxis, 'Xlim', px2, 'Ylim', py2);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 36
		px=get(gh.spineGUI.initialaxis, 'Xlim');
		px2=[px(1)+round(0.9*(px(2)-px(1))/2) px(2)];
		py=get(gh.spineGUI.initialaxis, 'Ylim');
		py2=[py(1)+round(0.9*(py(2)-py(1))/2) py(2)];
		set(gh.spineGUI.initialaxis, 'Xlim', px2, 'Ylim', py2);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 49 % show quadrant 1
		set(gh.spineGUI.initialaxis, 'Xlim', [1 round(1.1*size(state.imageProc.spine.initialImage,2)/2)], ...
			'Ylim', [1 round(1.1*size(state.imageProc.spine.initialImage, 1)/2)]);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 50 % show quadrant 2
		set(gh.spineGUI.initialaxis, 'Xlim', [round(0.9*size(state.imageProc.spine.initialImage,2)/2) size(state.imageProc.spine.initialImage,2)], ...
			'Ylim', [1 round(1.1*size(state.imageProc.spine.initialImage, 1)/2)]);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 51 % show quadrant 3
		set(gh.spineGUI.initialaxis, 'Xlim', [1 round(1.1*size(state.imageProc.spine.initialImage,2)/2)], ...
			'Ylim', [round(0.9*size(state.imageProc.spine.initialImage,1)/2) size(state.imageProc.spine.initialImage,1)]);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
	case 52 % show quadrant 4
		set(gh.spineGUI.initialaxis, 'Xlim', [round(0.9*size(state.imageProc.spine.initialImage,2)/2) size(state.imageProc.spine.initialImage,2)], ...
			'Ylim', [round(0.9*size(state.imageProc.spine.initialImage,1)/2) size(state.imageProc.spine.initialImage,1)]);
		if state.imageProc.spine.autoTransferTop
			transferImages;
		end
		
	case 20 % ctrl + t Load new image to Top spine analysis
		try
			saveSpineDataToExcel;
		catch
		end
		loadSpineImage;
	case 2 % ctrl + b Load new image to Bottom spine analysis
		loadSpineImageBot;
	case 5 % ctrl + e select spines
		if state.imageProc.spine.manual
			selectSpinesManually(gh.spineGUI.mainAxes);
		else
		end
	case 101 % e select spines
		if state.imageProc.spine.manual
			selectSpinesManually(gh.spineGUI.mainAxes);
		else
		end
	case 102 % f get FWHM
		selectSpinesManually_HW;
	case 103 % do spine analysis based on line
		selectSpinesManually_line;
		
	case 113 % q goto to frame 1
		state.imageProc.spine.maxFlag=0;
		state.imageProc.spine.currentSpineFrame = 1;
		updateGUIByGlobal('state.imageProc.spine.currentSpineFrame');
		set(state.imageProc.spine.imagehandle, 'CData', state.imageProc.spine.initialImage(:,:,state.imageProc.spine.currentSpineFrame));
		state.imageProc.spine.ZCurrentTop = state.imageProc.spine.ZStartTop + (state.imageProc.spine.currentSpineFrame-1)*state.imageProc.spine.zStepSizeTop;
		updateGUIByGlobal('state.imageProc.spine.ZCurrentTop');
		
		if state.imageProc.spine.autoTransferTop
			if state.imageProc.spine.auto
				showBWImage;
				return
			end
			
			if state.imageProc.spine.eraseSpines
				try
					removeLines;
				end
			end
			state.imageProc.spine.topImage=1;
			updateGuiByGlobal('state.imageProc.spine.topImage');
		
			state.imageProc.spine.bottomImage=1;
			updateGuiByGlobal('state.imageProc.spine.bottomImage');
			transferImages
		end

	case 115 % s advance 1 frame
		state.imageProc.spine.maxFlag=0;		
		state.imageProc.spine.currentSpineFrame = min(round(state.imageProc.spine.currentSpineFrame)+1, size(state.imageProc.spine.initialImage,3));
		updateGUIByGlobal('state.imageProc.spine.currentSpineFrame');
		set(state.imageProc.spine.imagehandle, 'CData', state.imageProc.spine.initialImage(:,:,state.imageProc.spine.currentSpineFrame));
		state.imageProc.spine.ZCurrentTop = state.imageProc.spine.ZStartTop + (state.imageProc.spine.currentSpineFrame-1)*state.imageProc.spine.zStepSizeTop;
		updateGUIByGlobal('state.imageProc.spine.ZCurrentTop');
		
		if state.imageProc.spine.autoTransferTop
			if state.imageProc.spine.auto
				showBWImage;
				return
			end
			
			if state.imageProc.spine.eraseSpines
				try
					removeLines;
				end
			end
			state.imageProc.spine.topImage=1;
			updateGuiByGlobal('state.imageProc.spine.topImage');
		
			state.imageProc.spine.bottomImage=1;
			updateGuiByGlobal('state.imageProc.spine.bottomImage');
			transferImages
		end
	case 97 % a to reverse 1 frame
		state.imageProc.spine.maxFlag=0;		
		state.imageProc.spine.currentSpineFrame = max(round(state.imageProc.spine.currentSpineFrame)-1, 1);
		updateGUIByGlobal('state.imageProc.spine.currentSpineFrame');
		set(state.imageProc.spine.imagehandle, 'CData', state.imageProc.spine.initialImage(:,:,state.imageProc.spine.currentSpineFrame));
		state.imageProc.spine.ZCurrentTop = state.imageProc.spine.ZStartTop + (state.imageProc.spine.currentSpineFrame-1)*state.imageProc.spine.zStepSizeTop;
		updateGUIByGlobal('state.imageProc.spine.ZCurrentTop');
		
		if state.imageProc.spine.autoTransferTop
			if state.imageProc.spine.auto
				showBWImage;
				return
			end
			
			if state.imageProc.spine.eraseSpines
				try
					removeLines;
				end
			end
			state.imageProc.spine.topImage=1;
			updateGuiByGlobal('state.imageProc.spine.topImage');
		
			state.imageProc.spine.bottomImage=1;
			updateGuiByGlobal('state.imageProc.spine.bottomImage');
			transferImages
		end

	case 12 % ctrl + l dendrte lengthe
		if state.imageProc.spine.manual
			getDendriteLength(gh.spineGUI.mainAxes);
			saveSpineDataToExcel;
		else
		end
	case 108 % l dendrte length
		if state.imageProc.spine.manual
			getDendriteLength(gh.spineGUI.mainAxes);
			saveSpineDataToExcel;
		else
		end
	case 104 % h -- measure area in a polygon
		getArea(gh.spineGUI.mainAxes);
		saveSpineDataToExcel;
	case 16 % ctrl + P open preview image
		loadPreviewSpine;
	case 19 % ctrl+s save the image
		saveMainAxis;

	case 122 % z dim image
		state.imageProc.spine.highPixelValue = round(state.imageProc.spine.highPixelValue*0.9);
	 	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
		set([gh.spineGUI.mainAxes gh.spineGUI.initialaxis], 'Clim', [state.imageProc.spine.lowPixelValue ...
			state.imageProc.spine.highPixelValue]);
	case 120 % x dim image
		state.imageProc.spine.highPixelValue = round(state.imageProc.spine.highPixelValue/0.9);
	 	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
		set([gh.spineGUI.mainAxes gh.spineGUI.initialaxis], 'Clim', [state.imageProc.spine.lowPixelValue ...
			state.imageProc.spine.highPixelValue]);
	case 117 % undo selection
		undoSpineSelection;
	otherwise
	end
catch
%	lasterr
end
