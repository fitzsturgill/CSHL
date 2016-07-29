function timerProcess_Physiology
	global state gh
	
	% put acq # and acq times in waves
	global physAcqTrace physAcqTime0 physAcqTime1 physCellVm0 physCellIm0 physCellVm1 physCellIm1
	physAcqTrace(state.files.lastAcquisition)=state.files.lastAcquisition;
	physAcqTime0(state.files.lastAcquisition)=state.phys.cellParams.minInCell0;
	physAcqTime1(state.files.lastAcquisition)=state.phys.cellParams.minInCell1;
	physCellVm0(state.files.lastAcquisition)=state.phys.cellParams.vm0;
	physCellVm1(state.files.lastAcquisition)=state.phys.cellParams.vm1;
	physCellIm0(state.files.lastAcquisition)=state.phys.cellParams.im0;
	physCellIm1(state.files.lastAcquisition)=state.phys.cellParams.im1;

	if ~isempty(state.internal.excelChannel) & state.files.autoSave
		try
			r=['r' num2str(25 + state.files.lastAcquisition)];
			ddepoke(state.internal.excelChannel, [r 'c30:' r 'c31'], [state.files.lastAcquisition state.epoch]);
			ddepoke(state.internal.excelChannel, [r 'c32'], state.phys.pulses.pulseSetName);
			ddepoke(state.internal.excelChannel, [r 'c33:' r 'c50'], [...
					state.cycle.lastPulseUsed0 ...
					state.cycle.lastPulseUsed1 ...
					state.phys.settings.extraGain0 ...
					state.phys.settings.extraGain1...
					state.phys.cellParams.minInCell0 ...
					state.phys.settings.currentClamp0 ...
					state.phys.cellParams.vm0 ...
					state.phys.cellParams.im0 ...
					state.phys.cellParams.rm0 ...
					state.phys.cellParams.rs0 ...
					state.phys.cellParams.cm0 ...
					state.phys.cellParams.minInCell1 ...
					state.phys.settings.currentClamp1 ...
					state.phys.cellParams.vm1 ...
					state.phys.cellParams.im1 ...
					state.phys.cellParams.rm1 ...
					state.phys.cellParams.rs1 ...
					state.phys.cellParams.cm1 ...
				]);
		catch
			disp('timerProcess_Physiology : unable to link to excel');
		end
	end

	try
		if state.analysis.active
			runTraceAnalyzer(1);
		end
	catch
		disp(['processPhysData : ' lasterr]);
		disp('	when doing trace analysis');
	end
		
	% call post analysis functions
	executeAnalysisFunctions;
	
	setPhysStatusString('');

 	% are they keeping data or killing it?
 	if ~state.phys.settings.keepInMemory
		for counter=1:length(state.phys.internal.newWaves)
			kill(state.phys.internal.newWaves{counter});
		end
	end		

	if ~state.cycle.loopingStatus
		set(gh.physControls.startButton, 'String', 'START');
		set([gh.physControls.startButton gh.scope.start], 'Enable', 'on');
	else
		set(gh.scope.start, 'Enable', 'on');
	end		
	
    try
	if state.db.conn~=0
		for counter=1:length(state.db.phys.wave_ids)
			state.db.phys.wave_id=state.db.phys.wave_ids(counter);
			state.db.phys.channel=state.db.phys.channels(counter);
			eval(['state.db.phys.currentClamp=state.phys.settings.currentClamp' num2str(counter) ';']);
			eval(['state.db.phys.inputGain=state.phys.settings.inputGain' num2str(counter) ';']);
			addRecordByTable('PhysAcq');
		end
		state.db.phys.wave_ids=[];
	    state.db.phys.channels=[];
	end
    catch
    end