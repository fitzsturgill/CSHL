function initializePhysDaq

	global state
	
	state.phys.daq.outputDevice = analogoutput('nidaq', state.phys.daq.outputBoardIndex);
	set(state.phys.daq.outputDevice, 'TriggerType', 'HwDigital');		
	set(state.phys.daq.outputDevice, 'SampleRate', state.phys.settings.outputRate*1000);
	state.phys.settings.outputRate=get(state.phys.daq.outputDevice, 'SampleRate')/1000;
	updateGuiByGlobal('state.phys.settings.outputRate');
	
	if state.phys.daq.auxOutputBoardIndex>0
		state.phys.daq.auxOutputDevice = analogoutput('nidaq', state.phys.daq.auxOutputBoardIndex);
		set(state.phys.daq.auxOutputDevice, 'TriggerType', 'HwDigital');
        % FS MOD  
        set(state.phys.daq.auxOutputDevice, 'HwDigitalTriggerSource', 'PFI0'); % PCI-6259 board is set as auxOutputDevice
        % END MOD
		set(state.phys.daq.auxOutputDevice, 'SampleRate', state.phys.settings.outputRate*1000);
		if ~timerGetActiveStatus('Imaging')
			syncNIDAQBoards(state.phys.daq.outputDevice, state.phys.daq.auxOutputDevice);	% Will sync the 2 board clocks.
		end
	end
	
	state.phys.daq.inputDevice = analoginput('nidaq', state.phys.daq.inputBoardIndex);
	set(state.phys.daq.inputDevice, 'TriggerType', 'HwDigital');				
	set(state.phys.daq.inputDevice, 'SamplesAcquiredFcn', {'processPhysData'});
	set(state.phys.daq.inputDevice, 'SampleRate', state.phys.settings.inputRate*1000);
    %FS MOD    
	set(state.phys.daq.inputDevice, 'InputType', 'SingleEnded');
    % END MOD
	state.phys.settings.inputRate=get(state.phys.daq.inputDevice, 'SampleRate')/1000;
	updateGuiByGlobal('state.phys.settings.inputRate');
    %% 
    % FS MOD
    % initialize multichannel acquisition
    mcAcqInitializePhysDaq;
    mcAcqInitializeOlfactometer; % adds DIO lines for olfactometer control
    % END MOD
    %%
    
	if timerGetActiveStatus('Imaging')
		syncNIDAQBoards(state.daq.focusOutput, state.phys.daq.outputDevice);	% Will sync the 2 board clocks.
	end
