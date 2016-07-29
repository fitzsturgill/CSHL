function setupInputChannels
	global state
	
	% delete active channels 
	delete(get(state.daq.focusInput, 'Channel'));
	delete(get(state.daq.grabInput, 'Channel'));
	delete(get(state.daq.pmtOffsetInput, 'Channel'));
	
	if state.acq.dualLaserMode==1 % if the lasers are on simulataneously then nothing special
		sampleFactor=1;
	elseif state.acq.dualLaserMode==2
		sampleFactor=2;	% if they are alternating, then double the number of acqs before trigger the trigger function
	else
		disp('	setupInputChannels needs more for lasermodes');
	end
	% add channels that are wanted but not active
	for channelCounter=1:state.init.maximumNumberOfInputChannels
		if getfield(state.acq, ['acquiringChannel' num2str(channelCounter)]) 
			channel=addchannel(state.daq.focusInput, channelCounter-1);
			channel.InputRange = [-10 10];
			channel=addchannel(state.daq.grabInput, channelCounter-1);
			channel.InputRange = [-10 10];
			channel=addchannel(state.daq.pmtOffsetInput, channelCounter-1);
			channel.InputRange = [-10 10];
		end
	end
	channel=addchannel(state.daq.pmtOffsetInput, 3);
	channel.InputRange = [-10 10];
	
	selectNumberOfStripes;	% select number of stripes based on # channels and resolution
	
	% GRAB acquisition: set up total acquisition duration
	actualInputRate = get(state.daq.grabInput, 'SampleRate');	
	state.internal.samplesPerLine = round(actualInputRate*state.acq.msPerLine);	
	state.internal.samplesPerFrame = state.internal.samplesPerLine*state.acq.linesPerFrame;

	% GRAB acquisition: set up action function trigger (1 per stripe)
	set(state.daq.grabInput, 'SamplesPerTrigger', sampleFactor*state.internal.samplesPerFrame*state.acq.numberOfFrames);
	set(state.daq.grabInput, 'SamplesAcquiredFcnCount', sampleFactor*state.internal.samplesPerFrame/state.internal.numberOfStripes); 

	% FOCUS acquisition: set up total acquisition duration
	actualInputRate = get(state.daq.focusInput, 'SampleRate');
	state.internal.samplesPerLineF = round(actualInputRate*state.acq.msPerLine);
	state.internal.samplesPerStripe = sampleFactor*state.internal.samplesPerLineF*state.acq.linesPerFrame/state.internal.numberOfStripes; 		 
	set(state.daq.focusInput, 'SamplesPerTrigger', ...
		state.internal.samplesPerStripe*state.internal.numberOfStripes*state.internal.numberOfFocusFrames);

	% FOCUS acquisition: set up action function trigger (1 per stripe)
	set(state.daq.focusInput, 'SamplesAcquiredFcnCount', state.internal.samplesPerStripe); 
	
	% PMT Offset: set up total acquisition duration 
	actualInputRate = get(state.daq.pmtOffsetInput, 'SampleRate');
	totalSamplesInputOffsets = 50*state.acq.samplesAcquiredPerLine;		% acquire 50 lines of Data
	set(state.daq.pmtOffsetInput, 'SamplesPerTrigger', totalSamplesInputOffsets);
