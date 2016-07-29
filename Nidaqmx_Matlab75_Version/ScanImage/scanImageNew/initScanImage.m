function initScanImage(userFile, analysisMode)

	global state gh
	
	if nargin<1
		userFile='';
		analysisMode=0;
	end
	
	if nargin==1
		if isnumeric(userFile) & ~isempty(userFile)
			analysisMode=userFile;
			userFile='';
		else
			analysisMode=0;
		end
	end
			
	if analysisMode
		h = waitbar(0, 'Starting ScanImage in Analysis Mode...', 'Name', 'ScanImage Analysis Initialization', 'WindowStyle', 'modal', 'Pointer', 'watch');
	else
		h = waitbar(0, 'Starting ScanImage...', 'Name', 'ScanImage Software Initialization', 'WindowStyle', 'modal', 'Pointer', 'watch');
	end
	gh.imageGUI = guihandles(imageGUI);
	gh.channelGUI = guihandles(channelGUI);
	gh.basicConfigurationGUI=guihandles(basicConfigurationGUI);
	gh.motorGUI =guihandles(motorGUI);
	gh.standardModeGUI=guihandles(standardModeGUI);
	gh.pcellControl = guihandles(pcellControl);
	gh.fieldAdjustGUI = guihandles(fieldAdjustGUI);

	set(gh.fieldAdjustGUI.scanRotationSlider, 'SliderStep', [5/360 15/360]);	% 5 degree changes for slider
	
	% Open the waitbar for loading
	
	waitbar(.1,h, 'Reading Initialization File...');
	openini('standard.ini');

	if analysisMode
		state.analysisMode=1;
		state.motor.motorOn=0;
	else
		state.analysisMode=0;
	end
	
	setStatusString('Initializing...');
	
	initPCellBoxSettingsManager;

	waitbar(.25,h, 'Creating Figures for Imaging');
	makeImageFigures;	% config independent...rleies only on the .ini file for maxNumberOfChannles.
	mp285Config;
	
	setStatusString('Initializing...');
	if ~analysisMode
		waitbar(.4,h, 'Setting Up Data Acquisition Devices...');
		
		setupImagingDaqs;
	end
	updateChannelFlags;
	if ~analysisMode
		updateDataForConfiguration;
	end
	updateKeepAllSlicesCheckMark; 
	updateCompostiteChannelSelections;
	
	initBlaster;
	
	if ~isempty(userFile)
		waitbar(.7,h, 'Reading User Settings...');
		openAndLoadUserSettings(userFile);
	else
		loadConfig;
	end
	
	makeConfigurationMenu;
	if ~analysisMode
		parkMirrors;
		closeShutter;
	end
	
	setStatusString('Initializing...');
	
	state.internal.status=0;
	applyChangesToOutput(1);

	waitbar(.9,h, 'Initialization Done');
	
	setStatusString('Ready to use');
	state.initializing=0;
	waitbar(1,h, 'Ready To Use');
	
	if analysisMode
		state.analysisMode=1;
		state.motor.motorOn=0;
 		set(get(gh.pcellControl.figure1, 'Children'), 'Enable', 'off');
 		set(get(gh.fieldAdjustGUI.figure1, 'Children'), 'Enable', 'off');

		set(get(gh.motorGUI.figure1, 'Children'), 'Enable', 'off');
		set(gh.standardModeGUI.focusButton, 'Enable', 'off')
		set(gh.standardModeGUI.grabOneButton, 'Enable', 'off')
        set(get(gh.basicConfigurationGUI.figure1, 'Children'), 'Enable', 'off')	
 %		set(get(gh.imageTracker.figure1, 'Children'), 'Enable', 'off')	
 %		set(gh.imageTracker.figure1, 'Visible', 'off')
 		set(gh.fieldAdjustGUI.figure1, 'Visible', 'off')
 		set(gh.pcellControl.figure1, 'Visible', 'off')
		set(gh.motorGUI.figure1, 'Visible', 'off');
		set(gh.blaster.figure1, 'Visible', 'off');
		set(state.internal.GraphFigure, 'Visible', 'off')
		set(state.internal.MaxFigure, 'Visible', 'off')
		set(state.internal.compositeFigure, 'Visible', 'off')
		set(gh.imageGUI.figure1, 'VIsible', 'off')		
 		initAvgAnalysis;
	else
		state.analysisMode=0;
	end

	close(h);