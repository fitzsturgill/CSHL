function getTriggerTime
% saves trigger time and calculates seconds since program startup

	global state

	state.internal.triggerTime = state.internal.dioTriggerTime;

	state.internal.triggerTimeString = clockToString(state.internal.triggerTime);
	state.internal.triggerTimeInSeconds = etime(state.internal.triggerTime, state.internal.startupTime);
	