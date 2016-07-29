function resetSlidersToMax
global state gh

set(gh.montageGUI.montageStartSlider, 'Min', 1 , 'Max', 1001, 'SliderStep', [.001 .001]);
set(gh.montageGUI.montageEndSlider, 'Min', 1 , 'Max', 1001, 'SliderStep', [.001 .001]);
set(gh.maxProjectionGUI.maxStartSlider, 'Min', 1 , 'Max', 1001, 'SliderStep', [.001 .001]);
set(gh.maxProjectionGUI.maxEndSlider, 'Min', 1 , 'Max', 1001, 'SliderStep', [.001 .001]);
set(gh.movieGUI.movieStartSlider, 'Min', 1 , 'Max', 1001, 'SliderStep', [.001 .001]);
set(gh.movieGUI.movieEndSlider, 'Min', 1 , 'Max', 1001, 'SliderStep', [.001 .001]);
set(gh.imageProcessingGUI.currentFrameSlider, 'Min', 1 , 'Max', 1001, 'SliderStep', [.001 .001]);
set(gh.imageProcessingGUI.totalFramesSlider, 'Min', 1 , 'Max', 1001, 'SliderStep', [.001 .001]);
	

