function getPower(level)
        global state
        state.pcell.pCellTestOutObj = analogoutput('nidaq',state.init.pockellBoardIndex);		
		state.pcell.pCellTestChannel = addchannel(state.pcell.pCellTestOutObj, state.init.pockellChannelIndex);	

		state.pcell.pCellTestInObj = analoginput('nidaq',state.init.acquisitionBoardIndex);
		state.pcell.pCellTestInChannel = addchannel(state.pcell.pCellTestInObj, 2);

		putsample(state.pcell.pCellTestOutObj, level);
        for counter=1:1
            pause(0.001);
			getsample(state.pcell.pCellTestInObj)
		end