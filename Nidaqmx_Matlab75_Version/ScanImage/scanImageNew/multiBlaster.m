function multiBlaster

global state;
global gh;

        setStatusString('Click on blaster locations');
		if state.internal.selectionChannel<=3
			figure(state.internal.GraphFigure(state.internal.selectionChannel));
		elseif state.internal.selectionChannel==4
			figure(state.internal.refFigure);
		end

    %locs=ginput;
    %state.blaster.multiLocs=locs;
    
    %for i=1:size(locs,1)
    %    addNewBlasterPosByScreenXY(locs(i, 1), locs(i,2));
    %end
    
    pos=ginput(1);
    
    while (length(pos)>0)
        state.blaster.multilocs(size(state.blaster.multilocs, 1)+1, 1)= pos(1);
        state.blaster.multilocs(size(state.blaster.multilocs,1), 2)= pos(2);
        addNewBlasterPosByScreenXY(pos(1), pos(2));
        updatereferenceimage
        pos=ginput(1);
    end
    
    setupBlasterConfigsByPos(22)
    makeBlasterConfigMenu
    updatereferenceimage
    
    %%%
   