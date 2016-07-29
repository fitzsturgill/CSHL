function computeTimeSeriesMeanV
% use Ctrl+s to compute
% Function will take the current image and compute the average intenisty over all the pixels in
% the frame.  The data is stored as a wave called imageProcMeanInt.

global state gh 

value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
frames = state.imageProc.cell.currentFrame{value}:state.imageProc.cell.numberOfFrames{value};
data = state.imageProc.cell.currentImage{value}(y:y1, x:x1, frames);

counter=1;
for i = 1:size(data,3)
	meanIntensity(counter,:) = mean(data(:,:,i));
	counter=counter+1;
end
% maxi=max(meanIntensity);
% meanIntensity=double(meanIntensity)/double(maxi);

c=wave('imageProcMeanInt', meanIntensity, 'note', ['Average for Image ' state.imageProc.cell.fileName{value}]);
if ~c
	global imageProcMeanInt
	imageProcMeanInt.data=meanIntensity;
	disp('imageProcMeanInt loaded with new data');
	if ~isempty(imageProcMeanInt.plot)
		parent=get(imageProcMeanInt.plot,'Parent');
		set(get(parent, 'XLabel'),'String', 'Time');
		set(get(parent, 'YLabel'),'String', 'Relative Intensity');
	end
else
	evalin('base', 'plot(imageProcMeanInt);');
	disp('imageProcMeanInt loaded with Intensity Data');
end
data = imageProcMeanInt.data;
name = state.imageProc.cell.fileName{value};

try
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column)] '''' ',name'''' );']);
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+1) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+length(data))] '''' ',data'''' );']);
	state.imageProc.spine.row = state.imageProc.spine.row+1;
	updateGUIByGlobal('state.imageProc.spine.row');
end


