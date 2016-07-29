function setupImagingDaqs
global state gh

% setups components of AO Objects that are config independent
    if state.analysisMode
        return
    end

% Mirror Data Output Acquisition (FOCUS)

	state.daq.focusOutput = analogoutput('nidaq',state.init.mirrorOutputBoardIndex);			
	addchannel(state.daq.focusOutput, state.init.XMirrorChannelIndex);
	addchannel(state.daq.focusOutput, state.init.YMirrorChannelIndex);
	set(state.daq.focusOutput, 'SampleRate', state.acq.outputRate);
	set(state.daq.focusOutput, 'TriggerType', 'HwDigital');					% Set to Trigger PFI6
	set(state.daq.focusOutput, 'RepeatOutput', state.internal.numberOfFocusFrames);

% Mirror Data Output Acquisition (GRAB)

	state.daq.grabOutput = analogoutput('nidaq',state.init.mirrorOutputBoardIndex);
	addchannel(state.daq.grabOutput, state.init.XMirrorChannelIndex);
	addchannel(state.daq.grabOutput, state.init.YMirrorChannelIndex);
	set(state.daq.grabOutput, 'SampleRate', state.acq.outputRate);
	set(state.daq.grabOutput, 'TriggerType', 'HwDigital');					% Set to Trigger PFI6

% Laser Parking Mirror Data Output Acquisition

	state.daq.parkMirrorsOutput = analogoutput('nidaq', state.init.mirrorOutputBoardIndex);			
	addchannel(state.daq.parkMirrorsOutput, state.init.XMirrorChannelIndex);
	addchannel(state.daq.parkMirrorsOutput, state.init.YMirrorChannelIndex);
	set(state.daq.parkMirrorsOutput, 'SampleRate', state.acq.outputRate);
	set(state.daq.parkMirrorsOutput, 'TriggerType', 'Immediate');					% Set to Trigger PFI6

	if state.pcell.pcellOn	% if using pockel cells
	% PCell output board for FOCUS
		state.daq.pcellFocusOutput = analogoutput('nidaq', state.pcell.pcellBoardIndex);			
		for counter=1:state.pcell.numberOfPcells
			addchannel(state.daq.pcellFocusOutput, counter-1); % getfield(state.pcell, ['pcellChannelIndex' num2str(counter)]));
		end

		for counter=1:state.pcell.numberOfPcells
			addchannel(state.daq.pcellFocusOutput, counter-1 + state.pcell.numberOfPcells); %getfield(state.pcell, ['pcellShutterIndex' num2str(counter)]));
		end

		set(state.daq.pcellFocusOutput, 'SampleRate', state.acq.outputRate);
		set(state.daq.pcellFocusOutput, 'TriggerType', 'HwDigital');					% Set to Trigger PFI6
		set(state.daq.pcellFocusOutput, 'RepeatOutput', state.internal.numberOfFocusFrames);
		
	% PCell output board for GRAB
		state.daq.pcellGrabOutput = analogoutput('nidaq', state.pcell.pcellBoardIndex);			
		for counter=1:state.pcell.numberOfPcells
			addchannel(state.daq.pcellGrabOutput, counter-1); %getfield(state.pcell, ['pcellChannelIndex' num2str(counter)]));
		end

		for counter=1:state.pcell.numberOfPcells
			addchannel(state.daq.pcellGrabOutput, counter-1 + state.pcell.numberOfPcells); %getfield(state.pcell, ['pcellShutterIndex' num2str(counter)]));
		end

		set(state.daq.pcellGrabOutput, 'SampleRate', state.acq.outputRate);
		set(state.daq.pcellGrabOutput, 'TriggerType', 'HwDigital');					% Set to Trigger PFI6
		syncNIDAQBoards(state.daq.focusOutput, state.daq.pcellFocusOutput);	% Will sync the 2 board clocks.
	end
	
	
%  Piezoelectrode focus control

	if state.piezo.usePiezo  %if using pz electrode
		state.piezo.Output=analogoutput('nidaq', state.piezo.pzBoardIndex);
		addchannel(state.piezo.Output, state.piezo.pzChannelIndex);
		set(state.piezo.Output, 'SampleRate', state.piezo.sampleRate);
	end

% % DIO Control of Shutter Opening 
	state.daq.shutterLine = addline(state.daq.dio, state.shutter.shutterLineIndex, 'out');
	closeShutter;

% INPUT OBJECTS BELOW
	
% GRAB input object
	state.daq.grabInput = analoginput('nidaq',state.init.acquisitionBoardIndex);
	set(state.daq.grabInput, 'TriggerType', 'HwDigital');										
	set(state.daq.grabInput, 'SampleRate', state.acq.inputRate);
	set(state.daq.grabInput, 'SamplesAcquiredFcn', {'makeFrameByStripes'});

% Focus Input object
	state.daq.focusInput = analoginput('nidaq',state.init.acquisitionBoardIndex);
	set(state.daq.focusInput, 'TriggerType', 'HwDigital');										% 6110E NI Board Set to Trigger PFI0
	set(state.daq.focusInput, 'SampleRate', state.acq.inputRate);
	set(state.daq.focusInput, 'SamplesAcquiredFcn', {'makeStripe'});

% PMT offsets
	state.daq.pmtOffsetInput = analoginput('nidaq', state.init.acquisitionBoardIndex);
	set(state.daq.pmtOffsetInput, 'TriggerType', 'Immediate');										% 6110E NI Board Set to Trigger PFI0
	set(state.daq.pmtOffsetInput, 'SampleRate', state.acq.inputRate);
	set(state.daq.pmtOffsetInput, 'SamplesAcquiredFcn', {});

