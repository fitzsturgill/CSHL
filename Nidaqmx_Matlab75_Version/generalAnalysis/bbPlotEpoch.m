function bbPlotEpoch(roi, epoch)
	
	if nargin==1
		global state
		epoch=state.epoch;
	end
	
	pulses=[1 4 5 8];
	colors={'blue', 'black', 'red', 'green'};
	
	for counter=1:length(pulses)
		pulse=pulses(counter);
		avgName=['e' num2str(epoch) 'p' num2str(pulse) 'c1r' num2str(roi) '_avg'];
		if ~iswave(avgName)
			wave(avgName, []);
		end
		if counter==1
			plot(avgName);
		else
			append(avgName);
		end
		setPlotProps(avgName, 'color', colors{counter});
	end
	
