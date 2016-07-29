function executeGrabOneCallback
	global state gh

	state.internal.looping=0;

	val=get(gh.standardModeGUI.grabOneButton, 'String');
		
	if strcmp(val, 'GRAB')
		if strcmp(get(gh.basicConfigurationGUI.figure1, 'Visible'), 'on') == 1
			beep;
			setStatusString('Close Configuration GUI');
			return
		end
		
		if ~savingInfoIsOK
			return
		end	

		timerSetActiveStatus('Imaging', 1);	

		timerCallPackageFunctions('FirstSetup', 'Imaging');
		timerCallPackageFunctions('Setup', 'Imaging');
		timerCallPackageFunctions('Start', 'Imaging');
		
		state.cycle.loopingStatus=0;

		if ~state.acq.externalTrigger
			diotrigger;
		end
	elseif strcmp(val, 'ABORT')
		timerCallPackageFunctions('Abort', 'Imaging');
	else
		disp('executeGrabOneCallback: Grab One button is in unknown state'); 	% BSMOD - error checking
	end
