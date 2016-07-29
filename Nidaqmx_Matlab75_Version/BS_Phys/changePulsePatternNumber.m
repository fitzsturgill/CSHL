function changePulsePatternNumber(n)

	global state

	if nargin<1
		n=state.phys.pulses.patternNumber;
	else
		state.phys.pulses.patternNumber=n;
		updateGUIByGlobal('state.phys.pulses.patternNumber');
	end

	fn=fieldnames(state.phys.pulses);

	for counter=1:length(fn)
		if findstr('List', fn{counter})
			fname=fn{counter}(1:end-4);
			if iscell(getfield(state.phys.pulses, fn{counter}))
				if length(getfield(state.phys.pulses, fn{counter}))>=n
					eval(['state.phys.pulses.' fname '=state.phys.pulses.' fn{counter} '{n};']);
				else
					eval(['state.phys.pulses.' fname '='''';']);
					eval(['state.phys.pulses.' fn{counter} '{n}='''';']);
				end
			else
				if length(getfield(state.phys.pulses, fn{counter}))>=n
					eval(['state.phys.pulses.' fname '=state.phys.pulses.' fn{counter} '(n);']);
				else
					if strcmp(fname, 'duration')
						eval(['state.phys.pulses.' fname '=1000;']);
						eval(['state.phys.pulses.'  fn{counter} '(n)=1000;']);
					else
						eval(['state.phys.pulses.' fname '=0;']);
						eval(['state.phys.pulses.'  fn{counter} '(n)=0;']);
					end
				end
			end			
			updateGUIByGLobal(['state.phys.pulses.' fname]);
		end
	end

	makePulsePattern(n);
