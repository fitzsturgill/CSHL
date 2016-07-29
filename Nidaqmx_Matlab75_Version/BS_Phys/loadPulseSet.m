function loadPulseSet(pname, fname)
	global state
	
	if nargin<2
		if ~isempty(state.phys.pulses.pulseSetPath)
			try
				cd(state.phys.pulses.pulseSetPath)
			catch
			end
		end
	
		[fname, pname]=uigetfile('*.mat', 'Choose pulse set');
	end

	if ~isnumeric(fname) & ~isempty(fname)
		if isempty(findstr(fname, '.'))
			fname=[fname '.mat'];
		end

		[fid, message]=fopen(fullfile(pname, fname));
		if fid<0
			disp(['loadPulseSet: Error opening pulseset file: ' message]);
			return
		end
		
		pulseSet=load(fullfile(pname, fname));
		pulseSet=pulseSet.pulseSet;
		fn=fieldnames(pulseSet);
				
		for counter=1:length(fn)
			if findstr('List', fn{counter})
				eval(['state.phys.pulses.' fn{counter} '= pulseSet.' fn{counter} ';']);
			end
		end

        if length(state.phys.pulses.patternRepeatsList)<length(state.phys.pulses.durationList)
            state.phys.pulses.patternRepeatsList(end+1:length(state.phys.pulses.durationList))=0;
            state.phys.pulses.patternISIList(end+1:length(state.phys.pulses.durationList))=0;
        end
		disp(['*** LOADED PULSE SET ' fullfile(pname, fname) ' ***']);
		
		setPhysStatusString('pulseSet loaded');
		state.phys.pulses.pulseSetChanged=0;

		periods=findstr(fname, '.');
		if any(periods)								
			fname=fname(1:periods(1)-1);
		end		
		state.phys.pulses.pulseSetName = fname;
		state.phys.pulses.pulseSetPath = pname;
		updateGUIByGlobal('state.phys.pulses.pulseSetName');		
		changePulsePatternNumber(1);
	else
		setPhysStatusString('Cannot open file');
	end
	
	