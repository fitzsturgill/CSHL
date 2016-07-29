function makeStripe(aiF, SamplesAcquired)
global state

% makeStripe.m*****
% Action Function
% Called during the focusMode.m script execution.
% Takes data from data acquisition engine and formats it into a proper intensity image.
%
% This function will take the datainput from the DAQ engine and remove the data for the
% lines that are acquired.  It will then bin the matrix along the columns to produce a final 1024 x 1024 image
%
% Written by: Thomas Pologruto
% Cold Spring Harbor Labs
% January 10, 2002


if state.internal.pauseAndRotate
	stopAndRestartFocus;
	return
end

try
	if state.internal.looping==1 & state.internal.stripeCounter==0;
		state.internal.secondsCounter=floor(state.internal.lastTimeDelay-etime(clock,state.internal.triggerTime));
		updateGUIByGlobal('state.internal.secondsCounter');
	end

	try
		stripeFinalData = double(getdata(state.daq.focusInput, state.internal.samplesPerStripe, 'native')); % Gets enoogh data for one stripe from the DAQ engine for all channels present
        %mean2(stripeFinalData) %TN 23jul05
        if state.internal.pauseAndRotate
			stopAndRestartFocus;
			return
		end
	catch
		return
	end

	inputChannelCounter = 0;

	startLine = 1 + state.acq.linesPerFrame/state.internal.numberOfStripes*state.internal.stripeCounter;
	stopLine = startLine + state.acq.linesPerFrame/state.internal.numberOfStripes - 1;

	for channelCounter = 1:state.init.maximumNumberOfInputChannels
		if state.internal.abortActionFunctions
			abortFocus;	
			return;
		end

		if state.acq.acquiringChannel(channelCounter)  % if statement only gets executed when there is a channel to focus.
			if getfield(state.acq, ['pmtOffsetAutoSubtractChannel' num2str(channelCounter)])
				offset=getfield(state.acq, ['pmtOffsetChannel' num2str(channelCounter)]); % get PMT offset for channel
			else
				offset=0;
			end
			
			inputChannelCounter = inputChannelCounter + 1;
			
			if state.acq.imagingChannel(channelCounter) | state.internal.composite
				if state.acq.dualLaserMode==1	% laser simultaneous
					tempStripe{channelCounter} = reshape(stripeFinalData(:, inputChannelCounter),  ...
						state.internal.samplesPerLineF, ...
						(state.acq.linesPerFrame/state.internal.numberOfStripes))' ...
						- offset;					% Extracts only Channel 1 Data
				elseif state.acq.dualLaserMode==2
					tempStripe{channelCounter} = reshape(stripeFinalData(:, inputChannelCounter),  ...
						state.internal.samplesPerLineF, ...
						(2*state.acq.linesPerFrame/state.internal.numberOfStripes))' ...
						- offset;					% Extracts only Channel 1 Data
				end

				if state.acq.bidi
					temp=tempStripe{channelCounter};
					temp(2:2:end,:)=fliplr(temp(2:2:end,:));
					temp = add2d(temp(:, state.internal.startDataColumnInLine:state.internal.endDataColumnInLine), state.acq.binFactor); %-offset;
					tempStripe{channelCounter}=temp;
				else		
					tempStripe{channelCounter} = add2d(tempStripe{channelCounter}(:, state.internal.startDataColumnInLine:state.internal.endDataColumnInLine), state.acq.binFactor); %-offset;
                end             
                
				if state.acq.dualLaserMode==1	% laser simultaneous
					state.acq.focusData{channelCounter}(startLine:stopLine,:) ...
						= tempStripe{channelCounter};
                         
					if state.acq.imagingChannel(channelCounter)
						% Displays the current images on the screen as they are acquired.
						set(state.internal.imagehandle(channelCounter), 'EraseMode', 'none', 'CData', tempStripe{channelCounter}, ...
							'YData', [startLine stopLine]); 		
					end					
					if state.internal.composite & state.internal.compositeChannelSelections(channelCounter)
						low=getfield(state.internal, ['lowPixelValue' num2str(channelCounter)]);
						hi=getfield(state.internal, ['highPixelValue' num2str(channelCounter)]);
						state.acq.compositeData(startLine:stopLine,:,state.internal.compositeChannelSelections(channelCounter))=...
							min(max(...
							(tempStripe{channelCounter} - low) / max(hi-low,1)...
							,0)...
							,1);
					end
				else
					if state.acq.dualLaserMode==2 %alternate by line
						numberOfLines=stopLine-startLine+1;
						state.acq.focusData{channelCounter}(startLine:stopLine,:) ...
							= tempStripe{channelCounter}(1:2:numberOfLines*2, :);
						state.acq.focusData{channelCounter+10}(startLine:stopLine,:) ...
							= tempStripe{channelCounter}(2:2:numberOfLines*2, :);
					end
					if state.acq.imagingChannel(channelCounter)
						% Displays the current images on the screen as they are acquired.
						set(state.internal.imagehandle(channelCounter), 'EraseMode', 'none', 'CData', tempStripe{channelCounter}(1:2:numberOfLines*2, :), ...
							'YData', [startLine stopLine]); 		
						set(state.internal.imagehandle(channelCounter+10), 'EraseMode', 'none', 'CData', tempStripe{channelCounter}(2:2:numberOfLines*2, :), ...
							'YData', [startLine stopLine]); 		
					end					
					if state.internal.composite 
						if state.internal.compositeChannelSelections(channelCounter)
							low=getfield(state.internal, ['lowPixelValue' num2str(channelCounter)]);
							hi=getfield(state.internal, ['highPixelValue' num2str(channelCounter)]);
							state.acq.compositeData(startLine:stopLine,:,state.internal.compositeChannelSelections(channelCounter))=...
								min(max(...
								(tempStripe{channelCounter}(1:2:numberOfLines*2, :) - low) / max(hi-low,1)...
								,0)...
								,1);
						end
						if state.internal.compositeChannelSelections(channelCounter+10)
							low=getfield(state.internal, ['lowPixelValue' num2str(channelCounter+10)]);
							hi=getfield(state.internal, ['highPixelValue' num2str(channelCounter+10)]);
							state.acq.compositeData(startLine:stopLine,:,state.internal.compositeChannelSelections(channelCounter+10))=...
								min(max(...
								(tempStripe{channelCounter}(2:2:numberOfLines*2, :) - low) / max(hi-low,1)...
								,0)...
								,1);
						end
					end
				end

			end
		end
	end

	if state.internal.composite			% BSMO 7/17/2
		set(state.internal.compositeImagehandle, 'EraseMode', 'none', 'CData', state.acq.compositeData(startLine:stopLine,:,:), ...
			'YData', [startLine stopLine]);
	end
	
	if state.internal.abortActionFunctions
		abortFocus;	
		return;
	end

	drawnow;

	state.internal.stripeCounter = state.internal.stripeCounter + 1; % increments the stripecounter to ensure proper image displays
	
	if  state.internal.stripeCounter == state.internal.numberOfStripes			
		state.internal.stripeCounter = 0;
		state.internal.focusFrameCounter = state.internal.focusFrameCounter + 1;
	end
		
	if state.internal.abortActionFunctions
		abortFocus;	
		return
	end

	if state.internal.focusFrameCounter + 1 == state.internal.numberOfFocusFrames
		abortFocus;
		return
	end
	
catch
	if state.internal.status~=2
		return
	end
	disp(['makeStripe:' lasterr]);
	if state.internal.abortActionFunctions
		abortFocus;	
		return
	end
end



		