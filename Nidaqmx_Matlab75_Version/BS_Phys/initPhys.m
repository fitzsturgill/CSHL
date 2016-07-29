function initPhys(userFile, analysisMode)
	if nargin<2
		userFile='';
		analysisMode=0;
	end
	
	global state gh
	
	if ~analysisMode
		h=waitbar(0,'Initializing Physiology');
	else
		h=waitbar(0,'Initializing Physiology in Analysis Mode');
	end	
	
	gh.physControls=guihandles(physcontrols);
	gh.physSettings=guihandles(physsettings);
	gh.scope=guihandles(scope);
	gh.pulseMaker=guihandles(pulseMaker);
    
    % FS MOD
    if ~isempty(state.phys.mcAcq.mcInputBoardIndex)
        gh.mcAcquisition=guihandles(mcAcquisition);
    end
    % end MOD
	
	
	waitbar(.1,h);
	
	openini('phys.ini');

    % FS MOD
    if ~isempty(state.phys.mcAcq.mcInputBoardIndex)
        mcAcqMakeProbeMenu;
        state.phys.mcAcq.displayData = repmat(NaN, 1, state.phys.mcAcq.totalChannels);
        state.phys.mcAcq.displayXData = NaN;
    end
    
    % end MOD
    
	state.phys.pulses.addCompList={''};
	state.phys.pulses.patternNameList={''};
	if analysisMode
		state.analysisMode=1;
	end
	
	waitbar(.2,h);
	
	waveo('physAcqTrace', repmat(nan, 1, 1000));
	waveo('physAcqTime0', repmat(nan, 1, 1000));
	waveo('physAcqTime1', repmat(nan, 1, 1000));	
	waveo('physCellRm0', repmat(nan, 1, 1000));
	waveo('physCellRs0', repmat(nan, 1, 1000));
	waveo('physCellCm0', repmat(nan, 1, 1000));
	waveo('physCellVm0', repmat(nan, 1, 1000));
	waveo('physCellIm0', repmat(nan, 1, 1000));
	waveo('physCellRm1', repmat(nan, 1, 1000));
	waveo('physCellRs1', repmat(nan, 1, 1000));
	waveo('physCellCm1', repmat(nan, 1, 1000));
	waveo('physCellVm1', repmat(nan, 1, 1000));
	waveo('physCellIm1', repmat(nan, 1, 1000));
	
	waitbar(.3,h);
	
	waveo('currentPulsePattern', 0);
	changePulsePatternNumber(1);
	global currentPulsePattern
	plot(currentPulsePattern);
	state.phys.internal.pulsePatternPlot=gcf;
	set(state.phys.internal.pulsePatternPlot, 'visible', 'off', ...
		'CloseRequestFcn', 'hideCurrentWindow', ...
		'Name', 'PULSE PATTERN', ...
		'NumberTitle', 'off', ...
		'Position', [400   656   321   144], ...
		'MenuBar', 'none');
	
	waitbar(.4,h);

	% initialize scope
	waveo('scopeOutput', 0, 'xscale', [0 1/state.phys.scope.outputRate]);
	waveo('scopeInput', 0, 'xscale', [0 1/state.phys.scope.inputRate]);
	waveo('scopeInputFit', []);
	global scopeInput scopeInputFit
	
	plot(scopeInput);
	state.phys.internal.scopeHandle=gcf;
	state.phys.internal.scopeAxisHandle=gca;
	state.phys.cellParams.breakInClock0=clock;
	state.phys.cellParams.breakInClock1=clock;	
	
	append(scopeInputFit);
	setPlotProps(scopeInputFit, 'color', 'red', 'LineWidth', 2);
	setPlotProps(scopeInput, 'LineWidth', 2);
	set(state.phys.internal.scopeHandle, 'Name', 'SCOPE', ...
		'NumberTitle', 'off', ...
		'CloseRequestFcn', 'hideScope', ...
		'MenuBar', 'none');
		
	waitbar(.5,h);

	changeChannelType([0 1]);
	initializeScopeDaq;
	makeScopeOutput;
	resetScopeDaq;
	waitbar(.6,h);

	if isfield(state, 'daq')
		if isfield(state.daq, 'focusOutput')		% there is a scanImage object created
			syncNIDAQBoards(state.daq.focusOutput, state.phys.daq.outputDevice);	% Will sync the 2 board clocks.
		end
	end
	
	if ~state.analysisMode
		state.phys.daq.triggerDevice = digitalio('nidaq', state.phys.daq.triggerBoardIndex);
		state.phys.daq.triggerLine = addline(state.phys.daq.triggerDevice, state.phys.daq.triggerLineIndex, 'out');
	
		initializePhysDaq;
		setupPhysDaqInputChannels;
	end
	
	if analysisMode
		set(gh.physControls.startButton, 'Enable', 'off')
 		set(get(gh.scope.figure1, 'Children'), 'Enable', 'off');
 		hideScope;
 		set(get(gh.pulseMaker.figure1, 'Children'), 'Enable', 'off');
	end
	waitbar(1,h);
	close(h);
