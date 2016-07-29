function abortCurrent
global state gh

% Function that checks the strings of FOCUS, GRAB, and LOOP and 
% aborts if they are running.
% This function is used in genericKeyPressFunction for the letter 'a'.
%
% Written By: Thomas Pologruto and Bernardo Sabatini
% Cold Spring Harbor Labs
% January 30, 2001

if state.internal.status==2
	abortFocus;
elseif state.internal.status==3
	abortGrab;
else
	if ~strcmp(get(gh.mainControls.startLoopButton, 'String'), 'LOOP')
		executeStartLoopCallback;
	elseif strcmp(get(gh.mainControls.focusButton, 'String'), 'ABORT')
		abortFocus;
	elseif strcmp(get(gh.mainControls.grabOneButton, 'String'), 'ABORT')
		abortGrab;
	end
end
	
