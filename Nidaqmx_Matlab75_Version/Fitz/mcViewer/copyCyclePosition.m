function copyCyclePosition(fromPosition, toPosition)
	global state

    
    if nargin < 2
        toPosition = length(state.cycle.delayList) + 1;
    end
    
    if nargin < 1
        fromPosition = state.cycle.displayCyclePosition;
    end

%     state.cycle.displayCyclePosition=position;
%     updateGuiByGlobal('state.cycle.displayCyclePosition');
        
        
	for counter=1:length(state.internal.cycleListNames)

        % copy cycle to new cycle position
        eval(['state.cycle.' state.internal.cycleListNames{counter} 'List(' num2str(toPosition) ') = state.cycle.'...
            state.internal.cycleListNames{counter} 'List(' num2str(fromPosition) ');']);
	end