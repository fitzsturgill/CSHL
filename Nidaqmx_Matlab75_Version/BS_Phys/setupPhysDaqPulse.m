function setupPhysDaqPulse

	global state

	if get(state.phys.daq.outputDevice, 'SamplesAvailable')>0
		flushPhysData;
	end

	if ~state.cycle.pulseToUse0 & ~state.cycle.pulseToUse1 %	no channels selected, use blank pulse on channel 0
		state.phys.internal.needNewChannels=1;
	end
	
	if state.phys.internal.needNewChannels
		delete(get(state.phys.daq.outputDevice, 'Channel'));
		added=0;
		
		if state.cycle.pulseToUse0
% 			addchannel(state.phys.daq.outputDevice, 0);
            ch = addchannel(state.phys.daq.outputDevice, 0); %FS MOD
            set(ch, 'OutputRange', [-10 10], 'UnitsRange', [-10 10]); %FS MOD
			added=1;
		end
		if state.cycle.pulseToUse1
% 			addchannel(state.phys.daq.outputDevice, 1);
            ch = addchannel(state.phys.daq.outputDevice, 1); %FS MOD
            set(ch, 'OutputRange', [-10 10], 'UnitsRange', [-10 10]); %FS MOD
			added=1;
		end
		
		if ~state.cycle.pulseToUse0 & ~state.cycle.pulseToUse1
% 			addchannel(state.phys.daq.outputDevice, 0);
            ch = addchannel(state.phys.daq.outputDevice, 0);    %FS MOD
            set(ch, 'OutputRange', [-10 10], 'UnitsRange', [-10 10]); %FS MOD            
			setPhysStatusString('NO OUTPUT SELECTED');
			disp('setupPhysDaqPulse : No ouput selected.  Turning on DA 0.');
		end
% 		if ~added
% 			error('setupPhysDaqPulse: No output channels selected');
% 		end
		state.phys.internal.needNewChannels=0;
		state.phys.internal.needNewOutputData=1;
	end
    
	if state.phys.internal.needNewOutputData
		output=[];
		for counter=0:1
			patternNum=eval(['state.cycle.da' num2str(counter) 'List(state.cycle.currentCyclePosition)']);
			RSPattern=[];
			if patternNum
				makePulsePattern(patternNum, 0);
				if getfield(state.phys.settings, ['channelType' num2str(counter)])>1
					if getfield(state.phys.settings, ['currentClamp' num2str(counter)])
						extraGain=getfield(state.phys.settings, ['extraGain' num2str(counter)]);
						gain=getfield(state.phys.settings, ['pAPerVOut' num2str(counter)]);
						if state.cycle.CCRCPulse
							makePulsePattern(state.cycle.CCRCPulse, 0);
							RSPattern=eval(['state.phys.pulses.pulsePattern' num2str(state.cycle.CCRCPulse)])';
						end
					else
						extraGain=getfield(state.phys.settings, ['extraGain' num2str(counter)]);
						gain=getfield(state.phys.settings, ['mVPerVOut' num2str(counter)]);
						if state.cycle.VCRCPulse
							makePulsePattern(state.cycle.VCRCPulse, 0);
							RSPattern=eval(['state.phys.pulses.pulsePattern' num2str(state.cycle.VCRCPulse)])';
						end
					end
				else
					extraGain=getfield(state.phys.settings, ['extraGain' num2str(counter)]);
					gain=1;
				end
				if isempty(output)
					output=(extraGain/gain)*eval(['state.phys.pulses.pulsePattern' num2str(patternNum)])';
				else
					pattern=eval(['state.phys.pulses.pulsePattern' num2str(patternNum)]);
					len=min(size(pattern, 2), size(output,1));
					output(1:len,end+1) = (extraGain/gain)*pattern(1:len)';
					if len<size(output,1)
						output(len+1:end,end)=output(len,end);		% fill with last point repeated
					end
				end
				if ~isempty(RSPattern)
					len=min(size(RSPattern,1),size(output,1));
					output(1:len)=output(1:len)+(RSPattern(1:len)/gain);
				end
			end
		end
		if ~state.cycle.pulseToUse0 & ~state.cycle.pulseToUse1
			disp('setupPhysDaqPulse : No ouput selected.  Using 1 second of blank output');
			output=zeros(round(1000*state.phys.settings.outputRate),1);
		end

		state.phys.daq.output=output;				% save the output data for next time
		state.phys.internal.needNewOutputData=0;
	end
	
	putdata(state.phys.daq.outputDevice, state.phys.daq.output);
	set(state.phys.daq.inputDevice, 'SamplesPerTrigger', round(size(state.phys.daq.output,1)*state.phys.settings.inputRate/state.phys.settings.outputRate));
	set(state.phys.daq.inputDevice, 'SamplesAcquiredFcnCount', round(size(state.phys.daq.output,1)*state.phys.settings.inputRate/state.phys.settings.outputRate));
    
    
    % FS MOD- copy Tom's modifications below
    % note the length of the multichannel acquisition is determined by the
    % phys output pulse... just like in scanimage
    if ~isempty(state.phys.mcAcq.mcInputBoardIndex) && state.cycle.mcOnList(state.cycle.currentCyclePosition)
        set(state.phys.mcAcq.mcInputDevice, 'SamplesPerTrigger', round(size(state.phys.daq.output,1)*state.phys.mcAcq.mcInputRate/state.phys.settings.outputRate));
        set(state.phys.mcAcq.mcInputDevice, 'SamplesAcquiredFcnCount', round(size(state.phys.daq.output,1)*state.phys.mcAcq.mcInputRate/state.phys.settings.outputRate));
    end
    
    
    
    % HANDLE AUX PHYS OUTPUT BOARD -- but only if the imaging board is not
    % doing it
    
    
 	if state.phys.daq.auxOutputBoardIndex 
    
        if state.phys.internal.needNewAuxOutputData & ~state.cycle.imageOnList(state.cycle.currentCyclePosition) %TN 09May05
    % TN 09May05
            chanNeeded=find(...
                [state.cycle.aux4List(state.cycle.currentCyclePosition) ...
                state.cycle.aux5List(state.cycle.currentCyclePosition) ...
                state.cycle.aux6List(state.cycle.currentCyclePosition) ...
                state.cycle.aux7List(state.cycle.currentCyclePosition)])+3;



			delete(get(state.phys.daq.auxOutputDevice, 'Channel'));
            
			if ~isempty(chanNeeded)
                
% 				addchannel(state.phys.daq.auxOutputDevice, chanNeeded);
% FS MOD
                addchannel(state.phys.daq.auxOutputDevice, chanNeeded - 4); % I want to access AO 0-3 for PCI-6259 board
% END		
			
				nPoints=size(state.phys.daq.output,1);
				state.phys.daq.auxOutput=zeros(nPoints, length(chanNeeded));
				
				counter=1;
				for channel=chanNeeded
					patternNum=eval(['state.cycle.aux' num2str(channel) 'List(state.cycle.currentCyclePosition);']);
					makePulsePattern(patternNum, 0);
					pattern=eval(['state.phys.pulses.pulsePattern' num2str(patternNum)]);
					pSize=size(pattern, 2);
					if nPoints > pSize
						pattern=[pattern repmat(pattern(end), 1, nPoints-pSize)];
                    end
					state.phys.daq.auxOutput(1:nPoints, counter)=pattern';
					counter=counter+1;
                end
                
				putdata(state.phys.daq.auxOutputDevice, state.phys.daq.auxOutput);
            end
            
%             if state.phys.mcAcq.olfEnabled && state.phys.mcAcq.olfMFCControl % FS MOD  
%                 addchannel(state.phys.daq.auxOutputDevice, chanNeeded - 4); % I want to access AO 0-3 for PCI-6259 board
%             end % FS MOD                
                
		end
	end
	state.phys.internal.outputNeedsUpdate=0;
	
			
