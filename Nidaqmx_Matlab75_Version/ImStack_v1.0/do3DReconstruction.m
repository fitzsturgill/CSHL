function do3DReconstruction
global gh state

f = figure('Position',  [505   415   543   537]);
state.imageProc.threedaxis=axes('YTickLabelMode', 'manual', 'XTickLabelMode', 'manual', 'ZTickLabelMode', 'manual', ...
	'YTickMode', 'manual', 'XTickMode', 'manual', 'ZTickMode', 'manual','Position', [0 0 1 1], 'Parent', f);


value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
	
	
a = state.imageProc.currentImage(y:y1,x:x1,state.imageProc.currentFrame:state.imageProc.numberOfFrames) > state.imageProc.highPixelValue;
v=double(a).*double(state.imageProc.currentImage(y:y1,x:x1,state.imageProc.currentFrame:state.imageProc.numberOfFrames));
r = size(v, 1);
c = size(v,2);
if state.imageProc.internal.autoShrink 
	v=genericBin(v,r/128,c/128,1);
end
v=uint16(v);
a=0;
p = patch(isosurface(v,2));
isonormals(v,p);
set(p,'FaceColor','blue','EdgeColor','none');
view(3); 
axis tight;
a=get(state.imageProc.threedaxis, 'DataAspectRatio');
state.imageProc.maxRatioX1 = a(1);
updateGUIByGlobal('state.imageProc.maxRatioX1');
state.imageProc.maxRatioY = a(2);
updateGUIByGlobal('state.imageProc.maxRatioY');
state.imageProc.maxRatioZ = a(3);
updateGUIByGlobal('state.imageProc.maxRatioZ');
