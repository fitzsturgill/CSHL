function timeSeriesROI
% This function comnputes and upodates the ROI sats for the
% roiAnalysis GUI in Imstack
global gh state
set(gh.roiAnalysis.figure1,'Pointer','watch');
value=get(gh.roiAnalysis.roiPopupmenu,'Value');
try
	counter=1;
	for i=state.imageProc.roiAnalysis.roiFrame:state.imageProc.roiAnalysis.roiEnd
		state.imageProc.roiAnalysis.imageData = state.imageProc.cell.currentImage{value}(:,:,i);
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
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			state.imageProc.roiAnalysis.timeSeriesData(counter)=mean(state.imageProc.roiAnalysis.roiData);
% 			state.imageProc.roiAnalysis.timeSeriesData(counter)=max(state.imageProc.roiAnalysis.roiData);
% 			state.imageProc.roiAnalysis.timeSeriesData(counter)=min(state.imageProc.roiAnalysis.roiData);
% 			state.imageProc.roiAnalysis.timeSeriesData(counter)=std(state.imageProc.roiAnalysis.roiData);
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			counter=counter+1;
		end
	end
	if evalin('base','exist(''TimeSeriesMean'')')
		if evalin('base','iswave(TimeSeriesMean);')
			evalin('base','TimeSeriesMean.data=state.imageProc.roiAnalysis.timeSeriesData;');
		else
			wave('TimeSeriesMean', state.imageProc.roiAnalysis.timeSeriesData);
		end
	else
		wave('TimeSeriesMean', state.imageProc.roiAnalysis.timeSeriesData);
	end
	disp('TimeSeriesMean loaded with Analysis Data as wave.');
end
set(gh.roiAnalysis.figure1,'Pointer','arrow');
try
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row+1) 'c' num2str(1) ':r' num2str(state.imageProc.spine.row+1) 'c' num2str(length(state.imageProc.roiAnalysis.timeSeriesData))] '''' ',state.imageProc.roiAnalysis.timeSeriesData'''' );']);
	state.imageProc.spine.row = state.imageProc.spine.row+1;
	updateGUIByGlobal('state.imageProc.spine.row');

end
