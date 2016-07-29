function mcAcqMakeProbeMenu
    global state gh
    
    children = get(gh.mcAcquisition.Probes, 'Children');
    
	if ~isempty(children)
		delete(children);
    end
    
    uimenu(gh.mcAcquisition.Probes, 'Label', 'Display Probe Order', 'Enable', 'on', 'Callback', 'disp(state.phys.mcAcq.mcChannelOrder)');
    
    state.phys.mcAcq.probeList = {'A1x16' 'A4x4' 'A2x2Tet' 'A1x16Poly2Std', 'A1x16Poly2s'};
    state.phys.mcAcq.channelOrderList  = {    [9 8 10 7 13 4 12 5 15 2 16 1 14 3 11 6 17 18] ... % 1x16
                                              [1 3 2 6 7 4 8 5 13 10 12 9 14 16 11 15 17 18] ... % 4x4
                                              [2 7 5 3 1 8 4 6 12 14 15 10 13 11 16 9 17 18]... %A2x2Tet
                                              [14 3 11 6 16 1 15 2 12 5 13 4 10 7 9 8 17 18]... %A1x16Poly2Std                                             
                                              [14 3 11 6 16 1 15 2 12 5 13 4 10 7 9 8 17 18]... %A1x16Poly2s                                              
                                              }; 
%         state.phys.mcAcq.channelOrderList  = {    [9 8 10 7 13 4 12 5 15 2 16 1 14 3 11 6 17 18] ... % 1x16
%                                               [1 3 2 6 7 4 8 5 13 10 12 9 14 16 11 15 17 18] ... % 4x4
%                                               [2 7 5 3 1 8 4 6 12 14 15 10 13 11 16 9 17 18]... %A2x2Tet
%                                               [3 14 6 11 1 16 2 15 5 12 4 13 7 10 8 9 17 18]... %A1x16Poly2Std                                             
%                                               [14 3 11 6 16 1 15 2 12 5 13 4 10 7 9 8 17 18]... %A1x16Poly2s                                              
%                                               }; 
%     
    
                                          
    for counter=1:length(state.phys.mcAcq.probeList)
        if counter==1
            uimenu(gh.mcAcquisition.Probes, 'Label', state.phys.mcAcq.probeList{1, counter}, 'Callback', 'mcAcqSelectProbeFromMenu' ...
                , 'Separator', 'on', 'Checked', 'on');
            state.phys.mcAcq.mcChannelOrder=state.phys.mcAcq.channelOrderList{1, 1};
            updateGuiByGlobal('state.phys.mcAcq.mcChannelOrder'); % no GUI yet I think...
        else
            uimenu(gh.mcAcquisition.Probes, 'Label', state.phys.mcAcq.probeList{1, counter}, 'Callback', 'mcAcqSelectProbeFromMenu');
        end
    end