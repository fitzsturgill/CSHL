function setAcquisitionParameters
global state gh

% setAcquisitionParameter.m****
% Function that sets the samplesAcquiredPerLine and inputRate when the user sets the 
% Fillfractiona dn the msPerLine.

state.acq.inputRate = 1250000;
updateGUIByGlobal('state.acq.inputRate');
		
switch state.acq.msPerLineGUI % 1 = 1 ms, 2 = 2ms, 3 = 4 ms, 4 = 8 ms
case 1
	state.acq.samplesAcquiredPerLine = 1024;
	updateGUIByGlobal('state.acq.samplesAcquiredPerLine');
	state.acq.outputRate = 125000;
	updateGUIByGlobal('state.acq.outputRate');
		
	switch state.acq.fillFractionGUI
	case 1 % fillFraction = 0.71234782608696
		state.acq.fillFraction =0.71234782608696;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .0011500;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 2 % fillFraction =  0.74472727272727
		state.acq.fillFraction = 0.7447272727272720000;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00110;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 3 % fillFraction = 0.78019047619048
		state.acq.fillFraction = 0.7801904761904800000;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .0010500;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 4 % fillFraction = 0.81920000000000
		state.acq.fillFraction = 0.8192000000000000000;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .0010;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 5 % fillFraction = 0.81920000000000
		state.acq.fillFraction = 0.8192000000000000000;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.fillFractionGUI = 4;
		updateGUIByGlobal('state.acq.fillFractionGUI');
		state.acq.msPerLine = .0010;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('Fill Fraction = .8192');
	case 6 % fillFraction = 0.81920000000000
		state.acq.fillFraction = 0.8192000000000000000;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.fillFractionGUI = 4;
		updateGUIByGlobal('state.acq.fillFractionGUI');
		state.acq.msPerLine = .0010;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('Fill Fraction = .8192');
	case 7 % fillFraction = 0.81920000000000
		state.acq.fillFraction = 0.8192000000000000000;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.fillFractionGUI = 4;
		updateGUIByGlobal('state.acq.fillFractionGUI');
		state.acq.msPerLine = .0010;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('Fill Fraction = .8192');
		
	otherwise 
	end
	
case 2
	state.acq.samplesAcquiredPerLine = 2048;
	updateGUIByGlobal('state.acq.samplesAcquiredPerLine');
	state.acq.outputRate = 125000;
	updateGUIByGlobal('state.acq.outputRate');
	
	switch state.acq.fillFractionGUI
	case 1 % fillFraction = 0.71234782608696
		state.acq.fillFraction = 0.71234782608696;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00230;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 2 % fillFraction =  0.74472727272727
		state.acq.fillFraction = 0.74472727272727;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00220;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 3 % fillFraction = 0.78019047619048
		state.acq.fillFraction = 0.78019047619048;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00210;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 4 % fillFraction = 0.81920000000000
		state.acq.fillFraction = 0.81920000000000;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .0020;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 5 % fillFraction = 0.86231578947368
		state.acq.fillFraction = 0.86231578947368;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00190;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 6 % fillFraction = 0.91022222222222
		state.acq.fillFraction = 0.91022222222222;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00180;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('');
	case 7  % fillFraction = 0.91022222222222
		state.acq.fillFraction = 0.91022222222222;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.fillFractionGUI = 6;
		updateGUIByGlobal('state.acq.fillFractionGUI');
		state.acq.msPerLine = .00180;
		updateGUIByGlobal('state.acq.msPerLine');
		setStatusString('Fill Fraction = .9102');
		
	otherwise 
	end
	
case 3
	
	setStatusString('');
	state.acq.samplesAcquiredPerLine = 4096;
	updateGUIByGlobal('state.acq.samplesAcquiredPerLine');
	state.acq.outputRate = 125000;
	updateGUIByGlobal('state.acq.outputRate');
		
	switch state.acq.fillFractionGUI
	case 1 % fillFraction = 0.71234782608696
		state.acq.fillFraction = 0.71234782608696;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00460;
		updateGUIByGlobal('state.acq.msPerLine');
	case 2 % fillFraction =  0.74472727272727
		state.acq.fillFraction = 0.74472727272727;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00440;
		updateGUIByGlobal('state.acq.msPerLine');
	case 3 % fillFraction = 0.78019047619048
		state.acq.fillFraction = 0.78019047619048;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00420;
		updateGUIByGlobal('state.acq.msPerLine');
	case 4 % fillFraction = 0.81920000000000
		state.acq.fillFraction = 0.81920000000000;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00400;
		updateGUIByGlobal('state.acq.msPerLine');
	case 5 % fillFraction = 0.86231578947368
		state.acq.fillFraction = 0.86231578947368;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00380;
		updateGUIByGlobal('state.acq.msPerLine');
	case 6 % fillFraction = 0.91022222222222
		state.acq.fillFraction = 0.91022222222222;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00360;
		updateGUIByGlobal('state.acq.msPerLine');
	case 7  % fillFraction = 0.96376470588235
		state.acq.fillFraction = 0.96376470588235;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .003400;
		updateGUIByGlobal('state.acq.msPerLine');
		
	otherwise 
	end
	
case 4
	setStatusString('');
	state.acq.samplesAcquiredPerLine = 8192;
	updateGUIByGlobal('state.acq.samplesAcquiredPerLine');
	state.acq.outputRate = 125000;
	updateGUIByGlobal('state.acq.outputRate');
	
	switch state.acq.fillFractionGUI
	case 1 % fillFraction = 0.71234782608696
		state.acq.fillFraction = 0.71234782608696;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00920;
		updateGUIByGlobal('state.acq.msPerLine');
	case 2 % fillFraction =  0.74472727272727
		state.acq.fillFraction = 0.74472727272727;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00880;
		updateGUIByGlobal('state.acq.msPerLine');
	case 3 % fillFraction = 0.78019047619048
		state.acq.fillFraction = 0.78019047619048;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00840;
		updateGUIByGlobal('state.acq.msPerLine');
	case 4 % fillFraction = 0.81920000000000
		state.acq.fillFraction = 0.81920000000000;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .0080;
		updateGUIByGlobal('state.acq.msPerLine');
	case 5 % fillFraction = 0.86231578947368
		state.acq.fillFraction = 0.86231578947368;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00760;
		updateGUIByGlobal('state.acq.msPerLine');
	case 6 % fillFraction = 0.91022222222222
		state.acq.fillFraction = 0.91022222222222;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .00720;
		updateGUIByGlobal('state.acq.msPerLine');
	case 7  % fillFraction = 0.96376470588235
		state.acq.fillFraction = 0.96376470588235;
		updateGUIByGlobal('state.acq.fillFraction');
		state.acq.msPerLine = .006800;
		updateGUIByGlobal('state.acq.msPerLine');
		
	otherwise 
	end
	
	
otherwise
end

state.internal.startDataColumnInLine = ...
	round(0.001*(state.acq.lineDelay+state.acq.mirrorLag)/state.acq.msPerLine*state.internal.samplesPerLine);
state.internal.endDataColumnInLine = state.internal.startDataColumnInLine + (state.acq.samplesAcquiredPerLine-1);	



