function computeTimeSeriesMean
% use Ctrl+a to compute
% Function will take the current image and compute the average intenisty over all the pixels in
% the frame.  The data is stored as a wave called imageProcMeanInt

global state gh imageProcMeanInt

value = get(gh.imageProcessingGUI.fileName, 'Value');
[x,x1,y,y1] = getCurrentAxisLimits(state.imageProc.internal.axis{value});
frames = state.imageProc.cell.currentFrame{value}:state.imageProc.cell.numberOfFrames{value};
data = state.imageProc.cell.currentImage{value}(y:y1, x:x1, frames);

counter=1;
for i = 1:size(data,3)
    meanIntensity(counter) = mean2(data(:,:,i));
    counter=counter+1;
end
maxi=max(meanIntensity);
meanIntensity=double(meanIntensity)/double(maxi);
if ~evalin('base', 'exist(''imageProcMeanInt'');')% ' %
	c=0;
else
	c=evalin('base', 'iswave(imageProcMeanInt);');% ' %
end
if ~c
    wave('imageProcMeanInt', meanIntensity, 'note', ['Average for Image ' state.imageProc.cell.fileName{value}]);
    imageProcMeanInt.data=meanIntensity;
    if ~isempty(imageProcMeanInt.plot)
        parent=get(imageProcMeanInt.plot,'Parent');
        set(get(parent, 'XLabel'),'String', 'Time');
        set(get(parent, 'YLabel'),'String', 'Relative Intensity');
    end
else
    if ~isempty(evalin('base','imageProcMeanInt.plot'))
        imageProcMeanInt.data=meanIntensity;
    else
        imageProcMeanInt.data=meanIntensity;
        evalin('base', 'plot(imageProcMeanInt);');
    end
end
data = imageProcMeanInt.data;
name = state.imageProc.cell.fileName{value};

try
    eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column)] '''' ',name'''' );']);
    eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+1) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+length(data))] '''' ',data'''' );']);
    state.imageProc.spine.row = state.imageProc.spine.row+1;
    updateGUIByGlobal('state.imageProc.spine.row');
end


