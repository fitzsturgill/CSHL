function applyAdvancedCyclePosition(position)
	global state

    if nargin==1
		state.cycle.currentCyclePosition=position;
		updateGuiByGlobal('state.cycle.currentCyclePosition');
	end
	
	timerSetActiveStatus('Imaging', state.cycle.imageOnList(state.cycle.currentCyclePosition));	
	if state.cycle.imageOnList(state.cycle.currentCyclePosition) % imaging is on
		try
			setupImagingCyclePosition;
		catch
			disp(['applyAdvancedCyclePosition : ' lasterr]);
		end
	end
	
	timerSetActiveStatus('Physiology', state.cycle.physOnList(state.cycle.currentCyclePosition));
    
    % FS MOD
    % I'm relying on the current cycle display variables thus I need to update them...
    % ex. #1  state.cycle.da0List and state.cycle.da0 are both display
    % variables whereas state.cycle.pulseToUse0 is a internal variable
    
    % ex. #2 state.timer.mcOlfValveList is for display whereas
    % state.timer.mcOlfValve is both for display and internal
	if state.phys.mcAcq.mcInputBoardIndex
        updateCycleDisplay(state.cycle.currentCyclePosition);
    end
    % END MOD
    
    
    
	if state.cycle.physOnList(state.cycle.currentCyclePosition) % physiology is on
		try
			if (state.cycle.da0List(state.cycle.currentCyclePosition) ~= state.cycle.pulseToUse0) ...
					| (state.cycle.da1List(state.cycle.currentCyclePosition) ~= state.cycle.pulseToUse1)
				state.phys.internal.needNewOutputData=1;
			end
			state.cycle.pulseToUse0=state.cycle.da0List(state.cycle.currentCyclePosition);
			state.cycle.pulseToUse1=state.cycle.da1List(state.cycle.currentCyclePosition);
			if (state.cycle.lastPulseUsed0 ~= state.cycle.pulseToUse0) | (state.cycle.lastPulseUsed1 ~= state.cycle.pulseToUse1)
				state.phys.internal.needNewOutputData=1;
				updateGUIByGlobal('state.cycle.pulseToUse0');
				updateGUIByGlobal('state.cycle.pulseToUse1');
			end
			
			channels=get(get(state.phys.daq.outputDevice, 'Channel'), 'HwChannel');
			if ~isnumeric(channels)
				for counter=1:length(channels)
					cList(counter)=channels{counter};
				end
			else
				cList=channels;
			end
			
			if length(cList)~=length(find([state.cycle.pulseToUse0 state.cycle.pulseToUse1]))
				state.phys.internal.needNewChannels=1;
			else
				if state.cycle.pulseToUse0 & isempty(cList==0)
					state.phys.internal.needNewChannels=1;
				end
				if state.cycle.pulseToUse1 & isempty(cList==1)
					state.phys.internal.needNewChannels=1;
				end
			end
			
            if state.phys.daq.auxOutputBoardIndex
                try
                    if ~all(state.cycle.lastUsedAuxPulses == [state.cycle.aux4List(state.cycle.currentCyclePosition) ...
                            state.cycle.aux5List(state.cycle.currentCyclePosition) ...
                            state.cycle.aux6List(state.cycle.currentCyclePosition) ...
                            state.cycle.aux7List(state.cycle.currentCyclePosition)])
   		                state.phys.internal.needNewAuxOutputData=1;
                    end
                catch
            		state.phys.internal.needNewAuxOutputData=1;
                end
                
%                 if state.phys.mcAcq.olfEnabled && state.cycle.mcOlfOn && state.phys.mcAcq.olfMFCControl % FS MOD  
%                     % last two aux outputs (from PCI-6259 board) are used
%                     % to control mass flow controllers, flow of diluent
%                     % stream increases during odor shunt to counteract
%                     % decrease in airflow of odor stream
%             		state.phys.internal.needNewAuxOutputData=1;                    % FS MOD
%                 end % FS MOD
            end
            

		catch
			disp(['applyAdvancedCyclePosition : ' lasterr]);
		end
    else
		state.cycle.pulseToUse0=0;
		state.cycle.lastPulseUsed0=0;
        state.phys.internal.abort=0;
    end
    
	state.cycle.nextTimeDelay=state.cycle.delayList(state.cycle.currentCyclePosition);
	updateGuiByGlobal('state.cycle.nextTimeDelay');

	state.cycle.lastPositionUsed=state.cycle.currentCyclePosition;
	
	
