function updateClim(channel)
	if nargin<1
		channel=[1 2 3 4 11 12 13 14];
	end
	global state

	
	for counter=channel
		low = getfield(state.internal, ['lowPixelValue' num2str(counter)]);
		high = getfield(state.internal, ['highPixelValue' num2str(counter)]);
		if ishandle(state.internal.axis(counter))
			set([state.internal.axis(counter) state.internal.maxaxis(counter)], 'CLim', [low high]);
			if state.internal.status==0 || state.internal.status==4 	% we are just sitting, so redraw window
				if state.internal.lastTaskDone==2	% we just focused
					set(state.internal.imagehandle(counter), 'EraseMode', 'none', 'CData', state.acq.focusData{counter}, ...
						'YData', [1 state.acq.linesPerFrame]); 	
				elseif state.internal.lastTaskDone==3 % we just grabbed
					frame=min(size(state.acq.acquiredData{counter},3),state.internal.reviewFrame);
					set(state.internal.imagehandle(counter), 'EraseMode', 'none', 'CData', state.acq.acquiredData{counter}(:,:,frame), ...
						'YData', [1 state.acq.linesPerFrame]); 	
				end
			end
					
		end
	end

	if state.internal.status==0 || state.internal.status==4
		state.acq.compositeData = (zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine, 3)); 	% BSMOD 7/17/2
		if state.internal.composite			% BSMOD 7/17/2
			for counter=1:3
				channel=find(state.internal.compositeChannelSelections==counter);
				if ~isempty(channel)
					channel=channel(end);
					if state.acq.acquiringChannel(mod(channel,10))
						low = getfield(state.internal, ['lowPixelValue' num2str(channel)]);
						high = getfield(state.internal, ['highPixelValue' num2str(channel)]);
						if state.acq.dualLaserMode==2 || (state.acq.dualLaserMode==1 && channel<=4)
							if state.internal.lastTaskDone==2	% we just focused
								state.acq.compositeData(:,:,counter)=...
									min(max(...
									(state.acq.focusData{channel} - low) / ...
									max(high,1)...
									,0)...
									,1);
							elseif state.internal.lastTaskDone==3 % we just grabbed
								state.acq.compositeData(:,:,counter)=...
									min(max(...
									(state.acq.acquiredData{channel}(:,:,end) - low) / ...
									max(high,1)...
									,0)...
									,1);
							end
						end
					end
				end
			end
			set(state.internal.compositeImagehandle, 'EraseMode', 'none', 'CData', ...
				state.acq.compositeData, 'YData', [1 state.acq.linesPerFrame]);
		end
	end