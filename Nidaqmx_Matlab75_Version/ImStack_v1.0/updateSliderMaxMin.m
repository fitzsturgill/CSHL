function updateSliderMaxMin
global gh state

value = get(gh.montageGUI.fileName, 'Value');
old=state.imageProc.cell.numberOfFrames{value};
state.imageProc.cell.numberOfFrames{value} = size(state.imageProc.cell.currentImage{value},3);

if state.imageProc.cell.numberOfFrames{value} == 1
	set(gh.montageGUI.montageStartSlider, 'Min', 1 , 'Max', 1.001, 'SliderStep', [1 1]);
	set(gh.montageGUI.montageEndSlider, 'Min', 1 , 'Max', 1.001, 'SliderStep', [1 1]);
else
	set(gh.montageGUI.montageStartSlider, 'Min', 1 , 'Max',  state.imageProc.cell.numberOfFrames{value}, ....
		'SliderStep', [1/(state.imageProc.cell.numberOfFrames{value}-1) 1/(state.imageProc.cell.numberOfFrames{value}-1) ]);

	set(gh.montageGUI.montageEndSlider, 'Min', 1 , 'Max',  state.imageProc.cell.numberOfFrames{value}, ....
		'SliderStep', [1/(state.imageProc.cell.numberOfFrames{value}-1) 1/(state.imageProc.cell.numberOfFrames{value}-1) ]);
end
		
value = get(gh.maxProjectionGUI.fileName, 'Value');
if state.imageProc.cell.numberOfFrames{value} == 1
	set(gh.maxProjectionGUI.maxStartSlider, 'Min', 1 , 'Max', 1.001, 'SliderStep', [1 1]);
	set(gh.maxProjectionGUI.maxEndSlider, 'Min', 1 , 'Max', 1.001, 'SliderStep', [1 1]);
else
	set(gh.maxProjectionGUI.maxStartSlider, 'Min', 1 , 'Max',  state.imageProc.cell.numberOfFrames{value}, ....
		'SliderStep', [1/(state.imageProc.cell.numberOfFrames{value}-1) 1/(state.imageProc.cell.numberOfFrames{value}-1) ]);

	set(gh.maxProjectionGUI.maxEndSlider, 'Min', 1 , 'Max',  state.imageProc.cell.numberOfFrames{value}, ....
		'SliderStep', [1/(state.imageProc.cell.numberOfFrames{value}-1) 1/(state.imageProc.cell.numberOfFrames{value}-1) ]);
end

value = get(gh.movieGUI.fileName, 'Value');
if state.imageProc.cell.numberOfFrames{value} == 1
	set(gh.movieGUI.movieStartSlider, 'Min', 1 , 'Max', 1.001, 'SliderStep', [1 1]);
	set(gh.movieGUI.movieEndSlider, 'Min', 1 , 'Max', 1.001, 'SliderStep', [1 1]);
else
	set(gh.movieGUI.movieStartSlider  , 'Min', 1 , 'Max',  state.imageProc.cell.numberOfFrames{value}, ....
		'SliderStep', [1/(state.imageProc.cell.numberOfFrames{value}-1) 1/(state.imageProc.cell.numberOfFrames{value}-1) ]);

	set(gh.movieGUI.movieEndSlider , 'Min', 1 , 'Max',  state.imageProc.cell.numberOfFrames{value}, ....
		'SliderStep', [1/(state.imageProc.cell.numberOfFrames{value}-1) 1/(state.imageProc.cell.numberOfFrames{value}-1) ]);
end

value = get(gh.imageProcessingGUI.fileName, 'Value');
if state.imageProc.cell.numberOfFrames{value} == 1
	set(gh.imageProcessingGUI.currentFrameSlider, 'Min', 1 , 'Max', 1.001, 'SliderStep', [1 1]);
	set(gh.imageProcessingGUI.totalFramesSlider, 'Min', 1 , 'Max', 1.001, 'SliderStep', [1 1]);
else
	
	set(gh.imageProcessingGUI.currentFrameSlider , 'Min', 1 , 'Max',  state.imageProc.cell.numberOfFrames{value}, ....
		'SliderStep', [1/(state.imageProc.cell.numberOfFrames{value}-1) 1/(state.imageProc.cell.numberOfFrames{value}-1) ]);

	set(gh.imageProcessingGUI.totalFramesSlider, 'Min', 1 , 'Max',  state.imageProc.cell.numberOfFrames{value}, ....
		'SliderStep', [1/(state.imageProc.cell.numberOfFrames{value}-1) 1/(state.imageProc.cell.numberOfFrames{value}-1) ]);
end
state.imageProc.cell.numberOfFrames{value}=old;
% get the size of the image (total number of frames)
updateGUIByGlobalCell('state.imageProc.numberOfFrames', value);
