function [maxImage, maxIndex]=extractMaxFromStack(channel, xb, yb)
	global state
	
	if nargin==0
		channel=1;
	end
	
	if nargin==3
		[maxImage, maxIndex]=max(state.acq.acquiredData{channel}(yb(1):yb(2), xb(1):xb(2), :), [], 3);
	else
		[maxImage, maxIndex]=max(state.acq.acquiredData{channel}, [], 3);
	end		