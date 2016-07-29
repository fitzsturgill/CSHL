function out=mp285ReadAnswer
	out=[];
	global state
	if length(state.motor.serialPortHandle) == 0
		disp(['MP285Talk: MP285 not configured']);
		return;
	end

	time=clock;
	done=0;
	while ~done
		n=get(state.motor.serialPortHandle,'BytesAvailable');
		if  n > 0
			temp=fread(state.motor.serialPortHandle,n); 
			out=[out temp];
			if temp(end)==13;
				done=1;
			end
			time=clock;
		end
		if etime(clock,time)>2
			disp('mp285ReadAnswer: Time out: no data in 2 secs');
			done=1;
		end
	end
		
