function removeLastStat
global gh state

try
	a = size(state.imageProc.spineData.allDataList,1);
	b = size(state.imageProc.spineData.allDataList,2);
	if state.imageProc.spineData.counter == 1
		return
	end
	
	state.imageProc.spineData.counter = state.imageProc.spineData.counter -1;
	
	for i = 1:(a-1)
		for j = 1:b	
			q{i,j} = state.imageProc.spineData.allDataList{i,j};
		end
	end
	state.imageProc.spineData.allDataList = q;
end
