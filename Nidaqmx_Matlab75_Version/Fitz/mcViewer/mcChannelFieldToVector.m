function out = mcChannelFieldToVector(f)
    % modified for mcAcquisition
    % gathers numeric or logical field values (field = f) into a vector
    global state
    if ~isfield(state.phys.mcAcq.channel, f)
        disp('***Error in mcChannelFieldToVector***');
        disp(['Field parameter == ' f]);
        return
    end
    if ~(isnumeric(state.phys.mcAcq.channel(1).(f)) || islogical(state.phys.mcAcq.channel(1).(f)))
            disp('***Error in mcChannelFieldToVector***');
            disp('not numeric or logical');
            f
            return
    end
    out = cat(2, state.phys.mcAcq.channel(1, :).(f));
    return