function initTimer(analysisMode, packages)

	global state gh
	
	if nargin<2
		analysisMode=0;
		packages=0;
	end
	
	if nargin==1
		if isnumeric(userFile) & ~isempty(userFile)
			analysisMode=userFile;
			userFile='';
		else
			analysisMode=0;
		end
    end
    
			
	gh.timerMainControls = guihandles(timerMainControls);
	set(gh.timerMainControls.statusString, 'String', 'Opening GUIS...')
	gh.advancedCycleGui=guihandles(advancedCycleGUI);
	set(gh.timerMainControls.statusString, 'String', 'Read timer.ini ...')
	openini('timer.ini');
	set(gh.timerMainControls.statusString, 'String', 'Read machineSpecific.ini ...')
	openini('machineSpecific.ini');
	setStatusString('Opening packages...');
	makeTimerPackagesMenu;

	if analysisMode
		state.analysisMode=1;
	else
		state.analysisMode=0;
	end
		
    if ~analysisMode
    	state.daq.dio = digitalio('nidaq', state.init.triggerBoardIndex);
    	state.daq.triggerLine = addline(state.daq.dio, state.init.triggerLineIndex, 'out');
    end

	state.internal.startupTime=clock;
	state.internal.startupTimeString=clockToString(state.internal.startupTime);
	updateHeaderString('state.internal.startupTimeString');

	state.internal.cycleListNames={};
	allfields=fieldnames(state.cycle);
	for counter=1:length(allfields)
		if ~isempty(findstr(allfields{counter}, 'List'))
			tag=allfields{counter};
			state.internal.cycleListNames{end+1}=tag(1:end-4);
		end
	end
	
	initNotebooks;
	if packages==1	% Imaging only
		state.initializing=1;
		timerSetActiveStatus('Imaging', 1);	
		state.cycle.imageOn=1;
		updateGUIByGlobal('state.cycle.imageOn');
		state.cycle.imageOnList(1)=1;
	elseif packages==2 % Phys only
		state.initializing=1;
		timerSetActiveStatus('Physiology', 1);	
		state.cycle.physOn=1;
		updateGUIByGlobal('state.cycle.physOn');
		state.cycle.physOnList(1)=1;
		seeGUI('gh.advancedCycleGui.figure1');
		if state.phys.daq.auxOutputBoardIndex & ~analysisMode
			syncNIDAQBoards(state.phys.daq.outputDevice, state.phys.daq.auxOutputDevice);	% Will sync the 2 board clocks.
		end

	elseif packages==3 % Both
		timerSetActiveStatus('Imaging', 1);	
		state.cycle.imageOn=1;
		updateGUIByGlobal('state.cycle.imageOn');
		state.cycle.imageOnList(1)=1;
		state.initializing=1;
		timerSetActiveStatus('Physiology', 1);	
		state.cycle.physOn=1;
		updateGUIByGlobal('state.cycle.physOn');
		state.cycle.physOnList(1)=1;
		seeGUI('gh.advancedCycleGui.figure1');
	end
	
	initTraceAnalysis;
	loadUserSettingsPath;
	loadUserSettings;
	updateCycleDisplay(1);

	if ~isempty(state.figurePath)
		try
			loadFiguresFromPath(state.figurePath);
		catch
			disp('initTimer : *** Error loading figure set ***');
		end	
	end
	
	if state.analysisMode
		set(gh.timerMainControls.loop, 'Enable', 'off');
		set(gh.timerMainControls.doOne, 'Enable', 'off');
		state.files.autoSave=0;
		updateGUIByGlobal('state.files.autoSave');
	end
	
	
	waveo('timerAcqTime', []);
	setStatusString('Ready to Use');
    
    %try
    %    TNSetPath;
    %catch
    %end
