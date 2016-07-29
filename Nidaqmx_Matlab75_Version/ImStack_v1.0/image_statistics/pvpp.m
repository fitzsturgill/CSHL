function pv=pvpp(chan)
	global state
	m=mean2(state.acq.acquiredData{chan});
	s=std2(state.acq.acquiredData{chan});

	m0=0;
	eval(['m0=state.acq.binFactor*state.acq.pmtOffsetChannel' num2str(chan) ';']);
	pv=s*s/(m-m0);
	disp(['offset = ' num2str(m0)]);
	disp(['mean = ' num2str(m)]);
	disp(['var = ' num2str(s*s)]);
	disp(['photons per pixel = ' num2str((m-m0)/pv)]);
	disp(['pvpp on channel ' num2str(chan) ' = ' num2str(pv)]);

	return