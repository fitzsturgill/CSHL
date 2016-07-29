function diotrigger(updateTrigger)
	if nargin<1
		updateTrigger=1;
	end
	
	global state
	
	state.internal.dioTriggerTime = clock;
	
	% Acquisition Board Trigger
	putvalue(state.daq.triggerLine, 1);			% Places an 'on' signal on the line initially
	putvalue(state.daq.triggerLine, 0); 			% Digital Trigger: Places a go signal (1 to 0 transition; FallingEdge) to 
												% the line to trigger the
												% ao1, ao2, & ai.
	if updateTrigger
		getTriggerTime;
	end
	updateHeaderString('state.internal.triggerTimeString');
	updateHeaderString('state.internal.triggerTimeInSeconds');
	