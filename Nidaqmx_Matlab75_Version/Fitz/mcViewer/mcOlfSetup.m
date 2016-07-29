function mcOlfSetup
    global state
    
    if state.phys.mcAcq.olfEnabled && state.cycle.mcOlfOn
        if state.phys.mcAcq.olfShuntEnabled
            set(state.phys.mcAcq.olfDevice, 'TimerPeriod', state.cycle.mcOlfDelay - state.phys.mcAcq.olfShuntDuration);
        else
             set(state.phys.mcAcq.olfDevice, 'TimerPeriod', state.cycle.mcOlfDelay);
        end
    end
    
    
    
    
%     if state.phys.mcAcq.olfEnabled && state.cycle.mcOlfOn
%         set(state.phys.mcAcq.olfDevice, 'TimerPeriod', state.cycle.mcOlfDelay);
%     end
    