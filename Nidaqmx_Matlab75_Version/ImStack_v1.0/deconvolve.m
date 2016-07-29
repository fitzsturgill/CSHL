function deconvolve
global gh state
% Will deconvolve the current image with a PSF...

value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});

if isempty(state.imageProc.PSF)
	[fname,pname]=uigetfile('*.tif', 'Select PSF');
	if isnumeric(fname)
		return
	else
		state.imageProc.PSF=opentif([pname fname]);
	end
end

answer = INPUTDLG('Enter Number of Iterations to Deconvolve','Iterations for Deconvolution');
if isempty(answer)
	answer=1;
else
	answer=str2num(answer{1});
end
state.imageProc.deconvolvedImage = [];
state.imageProc.deconvolvedImage=deconvolve3D(state.imageProc.currentImage(y:y1,x:x1,:),state.imageProc.PSF,answer);
loadImageFromArray('state.imageProc.deconvolvedImage');
