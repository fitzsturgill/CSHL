function roiStats
% This function comnputes and upodates the ROI sats for the
% roiAnalysis GUI in Imstack
global gh state

% Do analysis
if ~isempty(state.imageProc.roiAnalysis.imageData) & ishandle(state.imageProc.roiAnalysis.patchHandle)
    BW = roipoly(state.imageProc.roiAnalysis.imageData,state.imageProc.roiAnalysis.XPatch,...
        state.imageProc.roiAnalysis.YPatch);
    classImage=class(state.imageProc.roiAnalysis.imageData);
    NewImage=double(BW).*double(state.imageProc.roiAnalysis.imageData);
    sizeImg=size(NewImage);
    pixels=sizeImg(1)*sizeImg(2);
    rowOfImg=reshape(NewImage,pixels,1);   % MAde a columns vector
    index=find(NewImage > 0);
    state.imageProc.roiAnalysis.roiData=NewImage(index);
    state.imageProc.roiAnalysis.roiSum=	sum(state.imageProc.roiAnalysis.roiData);
    state.imageProc.roiAnalysis.roiMean=mean(state.imageProc.roiAnalysis.roiData);	
    state.imageProc.roiAnalysis.roiMedian=median(state.imageProc.roiAnalysis.roiData);
    state.imageProc.roiAnalysis.roiMax=max(state.imageProc.roiAnalysis.roiData);	
    state.imageProc.roiAnalysis.roiMin=min(state.imageProc.roiAnalysis.roiData);
	state.imageProc.roiAnalysis.roiNumber=length(state.imageProc.roiAnalysis.roiData);
	state.imageProc.roiAnalysis.roiStd=std(state.imageProc.roiAnalysis.roiData);
    updateGUIByGlobal('state.imageProc.roiAnalysis.roiSum');
    updateGUIByGlobal('state.imageProc.roiAnalysis.roiMean');
    updateGUIByGlobal('state.imageProc.roiAnalysis.roiMedian');
    updateGUIByGlobal('state.imageProc.roiAnalysis.roiMax');
    updateGUIByGlobal('state.imageProc.roiAnalysis.roiMin');
	updateGUIByGlobal('state.imageProc.roiAnalysis.roiNumber');
	updateGUIByGlobal('state.imageProc.roiAnalysis.roiStd');
    %     figure;
    %     imagesc(NewImage); 
end
