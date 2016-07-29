function addNewBlasterPosByScreenXY(x, y)

global state;

   
        state.blaster.displayPos=length(state.blaster.indexList)+1;
        updateGuiByGlobal('state.blaster.displayPos');
        state.blaster.XList(state.blaster.displayPos)=0;
		state.blaster.YList(state.blaster.displayPos)=0;
        %updatePositionDisplay;
    

    try
        x=round(x);
		y=round(y);
		
		index = round(y-1) * state.internal.lengthOfXData + round(state.internal.lengthOfXData*(state.internal.fractionStart + x*state.internal.fractionPerPixel));
		
		state.blaster.X = state.acq.rotatedMirrorData(index,1);
		state.blaster.Y = state.acq.rotatedMirrorData(index,2);
		
		updateGuiByGlobal('state.blaster.Y');
		updateGuiByGlobal('state.blaster.X');
		
		state.blaster.XList(state.blaster.displayPos)=state.blaster.X;
		state.blaster.YList(state.blaster.displayPos)=state.blaster.Y;
		state.blaster.indexXList(state.blaster.displayPos)=x;
		state.blaster.indexYList(state.blaster.displayPos)=y;
		state.blaster.indexList(state.blaster.displayPos)=index;
		
		setStatusString('');
        %state.blaster.currentConfig
        %if state.blaster.active & any(state.blaster.displayPos==state.blaster.allConfigs{state.blaster.currentConfig, 2}(:,1))	% the position is used in the current config
		%	state.internal.needNewRepeatedMirrorOutput=1;
		%	applyChangesToOutput;
        %end
    catch
		setStatusString('ERROR');
		disp(['blaster selection : ' lasterr]);
    end
	
	try
		%updateReferenceImage
    catch
		%disp(['blaster selection (updateReferenceImage): ' lasterr]);
    end
	%set(hObject, 'Enable', 'on');