function mcAcqInitializeOlfactometer

    global state
   
    % add 4 lines to triggerDevice to control PC-16 valve controller, 4
    % lines allow control of all 16 valves
    % line 0 for the scanimage triggerDevice is already utilized as the
    % central trigger signal which initiates acquisition.
    % Therefore I'm using lines 1 - 4 (posssibly more in future)
    state.phys.mcAcq.olfDevice = digitalio('nidaq', state.phys.mcAcq.olfBoardIndex);
    set(state.phys.mcAcq.olfDevice, 'TimerFcn', @mcAcqOlf_callback);

    state.phys.mcAcq.olfLines = addline(state.phys.mcAcq.olfDevice, state.phys.mcAcq.olfHwLines, state.phys.mcAcq.olfPort, 'out');
    if state.phys.mcAcq.olfShuntEnabled
        state.phys.mcAcq.olfShuntLine = addline(state.phys.mcAcq.olfDevice, state.phys.mcAcq.olfShuntHwLine, state.phys.mcAcq.olfPort, 'out');
    end
    state.phys.mcAcq.olfValveStatus = 0;  % initialize this so that valvce turns on with first acquisition
    state.phys.mcAcq.olfValveCode = logical([...
        1 0 0 0;... ch1
        0 1 0 0;... ch2
        1 1 0 0;... ch3
        0 0 1 0;... ch4
        1 0 1 0;... ch5
        0 1 1 0;... ch6
        1 1 1 0;... ch7
        0 0 0 1;... ch8
        1 0 0 1;... ch9
        0 1 0 1;... ch10
        1 1 0 1;... ch11
        0 0 1 1;... ch12
        1 0 1 1;... ch13
        0 1 1 1;... ch14
        1 1 1 1;... ch15
        0 0 0 0;... none or 16 if auto memory mode is utilized        
        ]);
    
    mcAcqValveSwitch(state.phys.mcAcq.olfDefaultValve); % turn dummy odor on