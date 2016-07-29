function mcAcqValveSwitch(index, shunt)
    % specific to "CODE" mode 
    % switches on valve using mcOlfDevice which controls 4 TTL lines
    if nargin < 2
        shunt = 0;
    end
    global state
    if ~state.phys.mcAcq.olfShuntEnabled
        putvalue(state.phys.mcAcq.olfDevice, state.phys.mcAcq.olfValveCode(index, :));
    elseif shunt
        putvalue(state.phys.mcAcq.olfDevice, [state.phys.mcAcq.olfValveCode(index, :) 1]);        
    else
        putvalue(state.phys.mcAcq.olfDevice, [state.phys.mcAcq.olfValveCode(index, :) 0]);
    end