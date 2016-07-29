function updateBlasterConfigLineDisplay
	global state
	state.blaster.linePos=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 1);
	state.blaster.linePat=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 2);
	state.blaster.lineWidth=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 3);
	state.blaster.linePower1=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 4);
	state.blaster.linePower2=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 5);
    try
    	state.blaster.lineTilerActive=state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 6);
    catch
        disp(['**** Incorrect blasterline format; using default: tiler off']);
        state.blaster.lineTilerActive=0;
    end
    updateGuiByGlobal('state.blaster.linePos');
	updateGuiByGlobal('state.blaster.linePat');
	updateGuiByGlobal('state.blaster.lineWidth');
	updateGuiByGlobal('state.blaster.linePower1');
	updateGuiByGlobal('state.blaster.linePower2');
	updateGuiByGlobal('state.blaster.lineTilerActive');

