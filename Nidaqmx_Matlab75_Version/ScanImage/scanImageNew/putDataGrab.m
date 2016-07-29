function putDataGrab
    global state
    
    putdata(state.daq.grabOutput, state.acq.repeatedMirrorData);			% Queues Data to engine for Board 2 (Mirrors)
    pcellChannels=0:2*state.pcell.numberOfPcells-1;
    
    if length(state.cycle.physOnList)>=state.cycle.currentCyclePosition & state.cycle.physOnList(state.cycle.currentCyclePosition) % BSMOD 08012005 Changes so that Imaging and Phys can both use the aux board
                                                                                                                                   %TN 03Aug05
        chanNeeded=[pcellChannels ...
                    find(...
                        [state.cycle.aux4List(state.cycle.currentCyclePosition) ...
                         state.cycle.aux5List(state.cycle.currentCyclePosition) ...
                         state.cycle.aux6List(state.cycle.currentCyclePosition) ...
                         state.cycle.aux7List(state.cycle.currentCyclePosition)])+3];
        
        delete(get(state.daq.pcellGrabOutput, 'Channel'));
        
        if isempty(chanNeeded)
            return
        end
            
        addchannel(state.daq.pcellGrabOutput, chanNeeded);
        rate=get(state.daq.pcellGrabOutput, 'SampleRate')/1000;
        nPoints=size(state.acq.pcellRepeatedOutput, 1);
        
        state.phys.daq.auxOutput=zeros(nPoints, length(chanNeeded));
        
        counter=1;
        for channel=chanNeeded
            if any(channel==pcellChannels)
                state.phys.daq.auxOutput(1:nPoints, counter)=state.acq.pcellRepeatedOutput(:, counter);
            else
                patternNum=eval(['state.cycle.aux' num2str(channel) 'List(state.cycle.currentCyclePosition);']);
                makePulsePattern(patternNum, 0, rate);
                pattern=eval(['state.phys.pulses.pulsePattern' num2str(patternNum)]);
                pSize=size(pattern, 2);
                if nPoints > pSize
                    pattern=[pattern repmat(pattern(end), 1, nPoints-pSize)];
                elseif pSize>nPoints
                    pattern=pattern(1:nPoints);
                end
                state.phys.daq.auxOutput(1:nPoints, counter)=pattern';
            end
            counter=counter+1;
        end
        putdata(state.daq.pcellGrabOutput, state.phys.daq.auxOutput);
    else              
        if size(get(state.daq.pcellGrabOutput, 'Channel'),1)~=size(pcellChannels,1)
            delete(get(state.daq.pcellGrabOutput, 'Channel'));
            addchannel(state.daq.pcellGrabOutput, pcellChannels);
        end
        putdata(state.daq.pcellGrabOutput, state.acq.pcellRepeatedOutput);		% Queues Data to engine for board 1 (Pockell Cell)
    end
    