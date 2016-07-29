function endFocus
global state gh

% endAcquisition.m*****
% Function called at the end of the Focus that will park the laser, close the shutter,
% reset the counters (internal), reset the currentMode, and make the 
% Grab One, Focus, and Loop buttons visible.
%
% Written By: Thomas Pologruto
% Cold Spring Harbor Labs
% January 26, 2001
	

	abortFocus;
	


	
	