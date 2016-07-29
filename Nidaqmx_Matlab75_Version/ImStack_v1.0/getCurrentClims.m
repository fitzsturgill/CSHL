function getCurrentClims
global state gh

clim = get(gca, 'CLim');
state.imageProc.generic.highPixelValue = clim(1,2);
updateGUIByGlobal('state.imageProc.generic.highPixelValue');

set(gh.genericLUTGUI.lowPixelSlider, 'SliderStep', [1/clim(2) 1/clim(2)], 'Max', 3*clim(2));
set(gh.genericLUTGUI.highPixelSlider, 'SliderStep', [1/clim(2) 1/clim(2)], 'Max', 3*clim(2));

state.imageProc.generic.lowPixelValue = clim(1,1);
updateGUIByGlobal('state.imageProc.generic.lowPixelValue');
