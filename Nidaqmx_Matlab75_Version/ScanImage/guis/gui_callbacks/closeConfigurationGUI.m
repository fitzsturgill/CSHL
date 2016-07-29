function closeConfigurationGUI
	global state

% Closes the configuration GUI and rebuilds DAQ devices if necesary
%
% Written By: Thomas Pologruto
% Cold Spring Harbor Labs
% January 26, 2001
% Modified: Bernardo Sabatini
% January 16, 2001
% Only resets devices if a change has been made to the configuration.

	try
		setStatusString('');
	catch
	end

	global gh
%	setfocus(gh.basicConfigurationGUI.samplesPerLine);
	if state.internal.configurationChanged 
		recordWindowPositions;
		updateDataForConfiguration;
	end
	state.internal.configurationChanged=0;

	hideGUI('gh.basicConfigurationGUI.figure1');
