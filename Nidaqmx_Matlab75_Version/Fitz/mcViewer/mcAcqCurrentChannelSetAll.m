function mcAcqCurrentChannelSetAll(fields)
    global state
    
    if nargin == 0 || isempty(fields)
        fields = fieldnames(state.phys.mcAcq.channel);
    end
    
   
    if ~iscell(fields)
        fields = {fields};
    end
    
    for i = 1:length(fields)
        field = fields{i};
        if strcmp(field, 'Name') || strcmp(field, 'MUAInclude') % skip these fields
            continue
        end
        for j = 1:length(state.phys.mcAcq.channel)
            state.phys.mcAcq.channel(j).(field) = state.phys.mcAcq.(['currentChannel' field]);
        end
    end