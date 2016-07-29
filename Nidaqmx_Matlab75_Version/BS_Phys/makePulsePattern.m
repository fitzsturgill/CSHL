function makePulsePattern(number, update, rate)
    global state

    if nargin<1
        number=state.phys.pulses.patternNumber;
    end
    if nargin<2
        update=1;
    end
    if nargin<3
        rate=state.phys.settings.outputRate;
    end

    
    if (state.phys.pulses.delayList(number) + state.phys.pulses.pulseWidthList(number) ...
        +state.phys.pulses.isiList(number)*(state.phys.pulses.numPulsesList(number)-1))...
            > state.phys.pulses.durationList(number)
        beep;
        setPhysStatusString('Error making pulse');
        error(['makePulsePattern: pulse #' num2str(number) ' duration is too short']);
    end
    
    data=state.phys.pulses.offsetList(number)...
         * ones(1, round(state.phys.pulses.durationList(number)*rate));
    
    if state.phys.pulses.pulseWidthList(number)>0
        for patternCounter=1:max(state.phys.pulses.patternRepeatsList(number), 1)
            start=state.phys.pulses.delayList(number) + (patternCounter-1)*state.phys.pulses.patternISIList(number);
            
            for counter=1:state.phys.pulses.numPulsesList(number)
                data(round(1+start*rate) ...
                     : round((start+state.phys.pulses.pulseWidthList(number))*rate)) ...
                    = state.phys.pulses.amplitudeList(number) + state.phys.pulses.offsetList(number);
                start=start + state.phys.pulses.isiList(number);
            end
        end	
    end
    
    for counter=str2num(state.phys.pulses.addCompList{number})
        if counter
            makePulsePattern(counter, update);
            addData=getfield(state.phys.pulses, ['pulsePattern' num2str(counter)]);
            len=min(length(addData), length(data));
            data(1:len)=data(1:len)+addData(1:len);
        end
    end
    eval(['state.phys.pulses.pulsePattern' num2str(number) '= data;']);
    
    if number==state.phys.pulses.patternNumber & update
        setWave('currentPulsePattern', 'data', data, 'xscale', [0 1/rate]);
    end
    
    if ~state.initializing
        if any(number==[state.cycle.da0List(state.cycle.currentCyclePosition)...
                        state.cycle.da1List(state.cycle.currentCyclePosition)])
            state.phys.internal.needNewOutputData=1;
        end        
        if any(number==[state.cycle.aux4List(state.cycle.currentCyclePosition) ...
                        state.cycle.aux5List(state.cycle.currentCyclePosition) ...
                        state.cycle.aux6List(state.cycle.currentCyclePosition) ...
                        state.cycle.aux7List(state.cycle.currentCyclePosition)])
            state.phys.internal.needNewAuxOutputData=1;
        end
    else
        state.phys.internal.needNewOutputData=1;
        state.phys.internal.needNewAuxOutputData=1;
    end
    
    
    

