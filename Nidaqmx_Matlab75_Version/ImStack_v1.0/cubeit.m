function h=cubeit(in,r,c,last)
global state gh

h=figure('position', [505   415   543   537]);
state.imageProc.cubeaxis = axes('Parent', h);
space = .5;
if state.imageProc.internal.autoShrink 
	if r > 128 | c > 128
	in=genericBin(in,r/128,c/128,1);
		r = 128;
		c=128;
		space = .1;
	end
end
meshgrid(1:space:r,1:space:c,1:space:last);
slice(in,[1 r],[1 c],[1 last], 'nearest');
colormap(gray);
shading(state.imageProc.cubeaxis, 'INTERP');
axis tight;
set(state.imageProc.cubeaxis, 'Clim', [state.imageProc.lowPixelValue state.imageProc.highPixelValue]);
set(state.imageProc.cubeaxis, 'ZTickLabel', '')
set(state.imageProc.cubeaxis, 'XTickLabel', '')
set(state.imageProc.cubeaxis, 'YTickLabel', '')
set(state.imageProc.cubeaxis, 'YTickMode', 'manual')
set(state.imageProc.cubeaxis, 'YTick', [])
set(state.imageProc.cubeaxis, 'ZTick', [])
set(state.imageProc.cubeaxis, 'XTick', [])
set(state.imageProc.cubeaxis, 'DrawMode', 'fast')
set(state.imageProc.cubeaxis, 'Position', [0 0 1 1]);
a=get(state.imageProc.cubeaxis, 'DataAspectRatio');
state.imageProc.maxRatioX1 = a(1);
updateGUIByGlobal('state.imageProc.maxRatioX1');
state.imageProc.maxRatioY = a(2);
updateGUIByGlobal('state.imageProc.maxRatioY');
state.imageProc.maxRatioZ = a(3);
updateGUIByGlobal('state.imageProc.maxRatioZ');


