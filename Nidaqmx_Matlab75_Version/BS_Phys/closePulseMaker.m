function closePulseMaker

	global state gh
	try
		hideGUI('gh.pulseMaker.figure1');
		set(state.phys.internal.pulsePatternPlot, 'Visible', 'off');
	catch
		disp('closePulseMaker: MATLAB must be closing.  Killing pulseMaker');
		delete(gh.pulseMaker.figure1);
		delete(state.phys.internal.pulsePatternPlot);
	end