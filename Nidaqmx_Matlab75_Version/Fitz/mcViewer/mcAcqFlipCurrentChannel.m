 function mcAcqFlipCurrentChannel(currentChannel)
    global state
    
    if nargin < 1
        currentChannel = state.phys.mcAcq.currentChannel;
    end
    
    if currentChannel > state.phys.mcAcq.totalChannels
        currentChannel = state.phys.mcAcq.totalChannels;
    end    
    state.phys.mcAcq.currentChannel = currentChannel;
    updateGUIByGlobal('state.phys.mcAcq.currentChannel');    
    fields = fieldnames(state.phys.mcAcq.channel);
    for i = 1:length(fields)
        field = fields{i};
        state.phys.mcAcq.(['currentChannel' field]) = state.phys.mcAcq.channel(currentChannel).(field);
        updateGUIByGlobal(['state.phys.mcAcq.currentChannel' field]);
    end
   