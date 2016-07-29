function mcOlfTrigger
    global state
    
    if state.phys.mcAcq.olfEnabled
        start(state.phys.mcAcq.olfDevice);
    end