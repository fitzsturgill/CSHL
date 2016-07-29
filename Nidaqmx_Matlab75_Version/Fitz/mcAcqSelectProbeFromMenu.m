function mcAcqSelectProbeFromMenu


    global state gh
    
    
    h=gcbo;
	children=get(gh.mcAcquisition.Probes, 'Children');
	index=find(children==h);
    
    name=get(children(index), 'Label');
    
    for i = 1:length(children)
        if i ~= index
            set(children(i), 'Checked', 'off')
        else
            set(children(i), 'Checked', 'on')
        end
    end
    
    
    index = find(strcmp(name, state.phys.mcAcq.probeList));
    
    state.phys.mcAcq.mcChannelOrder = state.phys.mcAcq.channelOrderList{1, index};
    updateGuiByGlobal('state.phys.mcAcq.mcChannelOrder');
    mcAcqUpdateChannelNames;
    mcUpdateFigures;

    