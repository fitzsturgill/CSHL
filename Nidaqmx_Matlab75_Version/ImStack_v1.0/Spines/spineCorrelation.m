function [meanR,stdDR] = spineCorrelation(rArray,distanceArray,binSize)
% do the spine correlation from excel data
global state

if nargin < 3
	binSize=2;
end

sizeRArray=size(rArray);
sizeDistanceArray=size(distanceArray);

newR=reshape(rArray, sizeRArray(1)*sizeRArray(2),1);
newD=reshape(distanceArray, sizeDistanceArray(1)*sizeDistanceArray(2),1);
[sortedD,ind] = sortrows(newD);
sortedR = newR(ind);

meanR=[];
stdDR=[];

for i=0:max(sortedD)
	ind=find(sortedD>i&sortedD<i+binSize);
	meanR = [meanR  mean(sortedR(ind))];
	stdDR = [stdDR  std(sortedR(ind))];
end
stdDR=stdDR/sqrt((sizeRArray(1)*sizeRArray(2))/2-1);

c=linspace(min(newD),max(newD),10);
figure
plot(meanR,'ro','color', 'blue');
hold on;
plot(stdDR,'ro','color', 'red');
figure
[out,n]=hist(newD,c);
plot(n,out/2);

number=length(meanR);

try
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+number)] '''' ',meanR'''' );']);
	state.imageProc.spine.row = state.imageProc.spine.row+1;
	updateGUIByGlobal('state.imageProc.spine.row');
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column+number)] '''' ',stdDR'''' );']);
	state.imageProc.spine.row = state.imageProc.spine.row+1;
	updateGUIByGlobal('state.imageProc.spine.row');
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column) ':r' num2str(state.imageProc.spine.row) 'c' num2str(state.imageProc.spine.column)] '''' ',binSize'''' );']);
	state.imageProc.spine.row = state.imageProc.spine.row+1;
	updateGUIByGlobal('state.imageProc.spine.row');
catch
	disp(['Uncable to connect to Excel. Start Spine Analyssi to establish a connection']);
end
