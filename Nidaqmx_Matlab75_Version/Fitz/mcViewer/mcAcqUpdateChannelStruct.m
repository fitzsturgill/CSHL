 function mcAcqUpdateChannelStruct(currentChannel)
    global state
    fields = fieldnames(state.phys.mcAcq.channel);
    for i = 1:length(fields)
        field = fields{i};
        state.phys.mcAcq.channel(state.phys.mcAcq.currentChannel).(field) = ...
            state.phys.mcAcq.(['currentChannel' field]);
    end