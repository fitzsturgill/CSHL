function setupPhysDaqInputChannels
% adds appropriate channels to state.phys.daq.inputDevice
	global state
	if state.analysisMode
		return
	end
	
	if size(get(state.phys.daq.inputDevice, 'Channel'),1)>0
		delete(state.phys.daq.inputDevice.Channel);
	end

	state.phys.settings.acquiredChannels=[];
	for counter=0:7
		if getfield(state.phys.settings, ['acq' num2str(counter)])
			channel=addchannel(state.phys.daq.inputDevice, counter);
			channel.InputRange = [-10 10];
            % FS MOD kludge, for respiration, 11-12-12
            if counter == 2
                disp('Kludge in setupPhysDaqInputChannels for respiration');
                channel.InputRange = [-1 1]; % avoid digitization of signal
            end
            % END MOD
			state.phys.settings.acquiredChannels(end+1)=counter;
		end
	end
	
	if isempty(get(state.phys.daq.inputDevice, 'Channel'))
		error('setupPhysDaqInputChannels: No input channels selected');
	else
		flushdata(state.phys.daq.inputDevice);
	end
	

	
