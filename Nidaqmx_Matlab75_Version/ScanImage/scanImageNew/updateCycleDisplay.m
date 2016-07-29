function updateCycleDisplay(position)
	global state

	if nargin<1
		position=state.cycle.displayCyclePosition;
	else
		if position>length(state.cycle.delayList)
			position=length(state.cycle.delayList)+1;
			create=1;
		else
			create=0;
		end
		state.cycle.displayCyclePosition=position;
		updateGuiByGlobal('state.cycle.displayCyclePosition');
	end
	
	for counter=1:length(state.internal.cycleListNames)
		if create
			if strcmp(state.internal.cycleListNames{counter}, 'frames')
				state.cycle.framesList(state.cycle.displayCyclePosition)=1;
			elseif strcmp(state.internal.cycleListNames{counter}, 'delay')
				state.cycle.delayList(state.cycle.displayCyclePosition)=10;
			elseif strcmp(state.internal.cycleListNames{counter}, 'repeats')
				state.cycle.repeatsList(state.cycle.displayCyclePosition)=1;
			elseif strcmp(state.internal.cycleListNames{counter}, 'zStepSize')
				state.cycle.zStepSizeList(state.cycle.displayCyclePosition)=1;
			elseif strcmp(state.internal.cycleListNames{counter}, 'numberOfZSlices')
				state.cycle.numberOfZSlicesList(state.cycle.displayCyclePosition)=1;
			else
				eval(['state.cycle.' state.internal.cycleListNames{counter} 'List(state.cycle.displayCyclePosition) = 0;']);
			end
		end
		eval(['state.cycle.' state.internal.cycleListNames{counter} '=state.cycle.' state.internal.cycleListNames{counter} 'List(state.cycle.displayCyclePosition);']);
		updateGuiByGlobal(['state.cycle.' state.internal.cycleListNames{counter}]);
	end
			
