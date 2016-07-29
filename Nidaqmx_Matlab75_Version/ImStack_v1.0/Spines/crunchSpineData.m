function crunchSpineData
global gh state

%Which Files To Crunch?
for i = 1:size(state.imageProc.spineData.allDataList,1)
	names{i} = state.imageProc.spineData.allDataList{i,1};
end
[s,v] = listdlg('PromptString','Select the Data to be crunched:', 'OKString', 'OK',...
                'SelectionMode','multiple',...
                'ListString', names, 'Name', 'Select a File');
if isempty(s)
	return
end

names = names(s);

totalColumns = size(state.imageProc.spineData.allDataList,2);

counter = 0;
for i = 1:size(state.imageProc.spineData.allDataList,1)
	for j = 1:length(names)
		if strcmp(state.imageProc.spineData.allDataList{i,1}, names(j))
			counter = counter+1;
			index(counter) = i;
			for k = 1:totalColumns
				newData{counter,k} = state.imageProc.spineData.allDataList{i,k};
			end
		end
	end
end

totalD = size(newData,1);

for m = 1:totalD
	name{m} = newData{m,1};
	data{m} = newData{m,2}(1:end);
	lengthHist{m} = newData{m,3};
	volHist{m} = newData{m,4};
end
meanRows = [1 3 5 7 9 11 13 15];
errorRows = [2 4 6 8 10 12 14 16];

count = 1;
for n = 1:length(data)
	meanVal(count,:) = data{n}(meanRows);
	errVal(count,:) = data{n}(errorRows);
	neurons(count) = data{n}(17);
	count = count+1;
end

neurons = sum(neurons);
errVal = std(meanVal); 
meanVal = mean(meanVal);

data{1}(meanRows) = meanVal;
data{1}(errorRows) = errVal;
data{1}(17) = neurons;
data = data{1};
name = name(1);

newlengthHist(:,1) = lengthHist{1}(:,1);
length(lengthHist);
for n = 1:length(lengthHist)
	newlengthHist(:,n+1) = lengthHist{n}(:,2);
end
tempHist = sum(newlengthHist(:,2:end),2);
newLengthHist = [newlengthHist(:,1) tempHist];

newVolHist(:,1) = volHist{1}(:,1);
length(lengthHist);
for n = 1:length(lengthHist)
	newVolHist(:,n+1) = volHist{n}(:,2);
end
tempHist = sum(newVolHist(:,2:end),2);
newVolHist = [newVolHist(:,1) tempHist];

crunchedData{1} = name;
crunchedData{2} = data;
crunchedData{3} = newLengthHist;
crunchedData{4} = newVolHist;

a = min(index);
state.imageProc.spineData.allDataList{a,1} = name;
for k = 2:size(crunchedData,2)
	state.imageProc.spineData.allDataList{a,k} = crunchedData{k};
end
state.imageProc.spineData.allDataList{a,1} = name{1};

newCounter = 1;
totalToRemove = size(index,2)-1;
TotalBefore = size(state.imageProc.spineData.allDataList,1);
NewTotal = TotalBefore - totalToRemove;
for k = 1:TotalBefore
	if ~iselement(index(2:end),k)
 		for j = 1:size(state.imageProc.spineData.allDataList,2)
 			final{newCounter,j} = state.imageProc.spineData.allDataList{k,j};
		end
		newCounter = newCounter+1;
	end
end

state.imageProc.spineData.allDataList = final;





				
