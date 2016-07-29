function histAndNoiseAnalysis
% This function displays the pixel histogram alonsdie with the image after 
% the moise is subtracted away;
global state gh
value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
img=state.imageProc.cell.currentImage{value}(:,:,state.imageProc.cell.currentFrame{value});
[state.imageProc.hist,state.imageProc.n,meanImg,stdImg,SEM]=processImageDataYalin(img,state.imageProc.parsing.threshold,state.imageProc.parsing.meanAn,...
	state.imageProc.parsing.stdAn);
meanImg
stdImg