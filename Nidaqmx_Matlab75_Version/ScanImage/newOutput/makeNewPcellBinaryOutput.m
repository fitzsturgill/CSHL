function makeNewPcellBinaryOutput
	global state gh
	
	state.internal.lineDelay = 0.001*state.acq.lineDelay/state.acq.msPerLine; % calculate fractional line delay
	pOut = zeros(state.internal.lengthOfXData ,1);					% Fill with zero for the flyback
	if state.acq.bidi
		state.internal.flybackDecimal = 1-state.acq.fillFraction-state.internal.lineDelay;
		pStart =  1+round(state.internal.lengthOfXData * ...
			(0.001*(state.acq.lineDelay+state.acq.mirrorLag-state.pcell.pcellDelay)/state.acq.msPerLine+state.internal.flybackDecimal/2));
		pEnd =  1+round(state.internal.lengthOfXData * ...
			(0.001*(state.acq.lineDelay+state.acq.mirrorLag-state.pcell.pcellDelay)/state.acq.msPerLine+state.acq.fillFraction+state.internal.flybackDecimal/2));
	else
		pStart = round(state.internal.lengthOfXData * ...
			0.001*(state.acq.lineDelay+state.acq.mirrorLag-state.pcell.pcellDelay)/state.acq.msPerLine);
		pEnd =  round(state.internal.lengthOfXData * ...
			(0.001*(state.acq.lineDelay+state.acq.mirrorLag-state.pcell.pcellDelay)/state.acq.msPerLine+state.acq.fillFraction));
	end
	pOut(max(pStart,1):min(pEnd,state.internal.lengthOfXData))=1;			% Fill with 1 for scanning portion
	state.acq.pcellSingleLineBinary=pOut;
	
	if state.acq.dualLaserMode==1
		state.acq.pcellBinaryOutput = repmat(pOut, [state.acq.linesPerFrame 1]); 						% Final Pockell Data for one frame
	elseif state.acq.dualLaserMode==2
		state.acq.pcellBinaryOutput = repmat([pOut' zeros(1, length(pOut))]', [state.acq.linesPerFrame 1]); 						% Final Pockell Data for one frame
		state.acq.pcellBinaryOutputComp = repmat([zeros(1, length(pOut)) pOut']', [state.acq.linesPerFrame 1]); 						% Final Pockell Data for one frame
	end		
	
