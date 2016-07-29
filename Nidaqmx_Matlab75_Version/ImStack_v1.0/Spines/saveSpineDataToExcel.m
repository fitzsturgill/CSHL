function saveSpineDataToExcel
global state

peaks=[];

if state.imageProc.spine.dendriteLength == 1
	number = [state.imageProc.spine.numberOfSpines];
else
	number= [state.imageProc.spine.numberOfSpines state.imageProc.spine.dendriteLength state.imageProc.spine.spineDensity];
end

% if state.imageProc.spine.dendriteLength == 1
% 	calc=0;
% else
% 	calc=2;
% end

calc=size(number,2)-1;
data = state.imageProc.spine.spineLengths;
columns = size(data,2)-1;

if state.imageProc.spine.topImage
	name = state.imageProc.spine.loadedFileNameTop;
elseif state.imageProc.spine.bottomImage
	name = state.imageProc.spine.loadedFileNameBot;
else
	name =''
end

space = ' ';

try
	size = [valueFromHeaderString('state.internal.triggerTimeInSeconds', state.imageProc.spine.headerTop) valueFromHeaderString('state.acq.zoomFactor', state.imageProc.spine.headerTop)];
catch
	size = [1 0];
end

if isempty(size(1)) 
    size=[0 0];
end

eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column)] '''' ',name'''' );']);
eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+1) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+2)] '''' ',size'''' );']);
eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+3) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+3+calc)] '''' ',number'''' );']);
if calc == 0
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+4) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+4)] '''' ',space'''' );']);
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+5) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+5)] '''' ',space'''' );']);
	calc = calc+2;
end
eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+4+calc) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+4+calc)] '''' ',space'''' );']);
peaks=[];
if ~isempty(peaks)
	peakLength=length(peaks)-1;
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+5+calc) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+5+calc+peakLength)] '''' ',peaks'''' );']);
else
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+5+calc) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+5+calc)] '''' ',space'''' );']);
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+6+calc) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+6+calc)] '''' ',space'''' );']);
	peakLength=0;
end

if ~isempty(data)
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+peakLength+7+calc) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+peakLength+7+columns+calc)] '''' ',data'''' );']);
else
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+peakLength+7+calc) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+peakLength+7+columns+calc)] '''' ',space'''' );']);
end
for counter=1:10
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+peakLength+7+calc+length(data)+counter) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+peakLength+7+columns+calc+length(data)+counter)] '''' ',space'''' );']);
end
