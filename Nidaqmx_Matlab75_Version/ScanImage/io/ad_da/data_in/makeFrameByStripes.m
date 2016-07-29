function makeFrameByStripes(ai, SamplesAcquired)
global state

% makeFrame.m*****
% Data Acquisition (SamplesAcquired) Action Function
% Used with the contMode.m script to update frames on the screen after each frame.
% Takes data from data acquisition engine and formats it into a proper intensity image.
%
% This function will take the datainput from the DAQ engine and remove the data for the
% lines that are acquired.  It will then bin the matrix along the columns to produce a final image
% 
% The image will update every frame on the screen as data is recorded.
% The data is stored in the cell array state.acq.acquiredData{X}(:,:,frames)(X = 1,2,3...) , where X is the channel Acquired the frames are indexed in the third dimension.
% 
% This action function also handles averaging over frames.
% 
% Written by: Thomas Pologruto and Bernardo Sabatini
% Cold Spring Harbor Labs
% January 31, 2001

% Write complete header string  for only the first frame
	
	if state.internal.stripeCounter==0
		state.internal.frameCounter = state.internal.frameCounter + 1;	% Increments the frameCounter to ensure proper image storage and display
		updateGUIByGlobal('state.internal.frameCounter');	% Updates the frame Counter on the main controls GUI.

		if state.internal.looping==1
			state.internal.secondsCounter=floor(state.internal.lastTimeDelay-etime(clock,state.internal.triggerTime));
		else
			state.internal.secondsCounter=floor(etime(clock,state.internal.triggerTime));
		end
		updateGuiByGlobal('state.internal.secondsCounter');
	end
	
	try
		if state.internal.frameCounter == state.acq.numberOfFrames & state.internal.stripeCounter==state.internal.numberOfStripes-1
			closeShutter;
			if state.acq.numberOfZSlices > 1
				startMoveStackFocus; 	% start movement - focal plane down one step
			end
		end

		lps=state.acq.linesPerFrame/state.internal.numberOfStripes;
		startLine =1 + state.internal.stripeCounter*lps ;
		stopLine= startLine+lps-1;

		if state.internal.keepAllSlicesInMemory
			if state.acq.averaging
				position =state.internal.zSliceCounter + 1;
			else
				position =(state.internal.frameCounter + state.internal.zSliceCounter*state.acq.numberOfFrames);
			end
		else
			if state.acq.averaging
				position=1;
			else
				position = state.internal.frameCounter;
			end
		end
	catch
		disp(['Error in makeFromByStripes (1) : ' lasterr]);
	end

	try
		if state.acq.dualLaserMode==1 % simultaneous lasers
			frameFinalData = double(getdata(state.daq.grabInput, state.internal.samplesPerFrame/state.internal.numberOfStripes, 'native')); % Gets enoogh data for one frame from the DAQ engine for all channels present
		elseif state.acq.dualLaserMode==2 % alternate by line
			frameFinalData = double(getdata(state.daq.grabInput, 2*state.internal.samplesPerFrame/state.internal.numberOfStripes, 'native')); % Gets enoogh data for one frame from the DAQ engine for all channels present
		end
	catch
		if ~state.internal.abortActionFunctions & state.internal.status==3	% surpress errors if we are aborting or we have already aborted
			disp(['Error in makeFromByStripes (getdata) : ' lasterr]);
			timerCallPackageFunctions('Abort', 'Imaging');
			return
		else
			abortGrab;
			return
		end
	end
	
	try
		inputChannelCounter = 0;
		if state.internal.abortActionFunctions
			abortGrab;
			return
        end
	
		for channelCounter = 1:state.init.maximumNumberOfInputChannels
			if state.acq.acquiringChannel(channelCounter) % if statemetnt only gets executed when there is a channel to acquire.
				inputChannelCounter = inputChannelCounter + 1;
				if getfield(state.acq, ['pmtOffsetAutoSubtractChannel' num2str(channelCounter)])
					offset=getfield(state.acq, ['pmtOffsetChannel' num2str(channelCounter)]);
				else
					offset=0;
				end

				if state.acq.dualLaserMode==1	% laser simultaneous
					tempImage = reshape(frameFinalData(:,inputChannelCounter), state.internal.samplesPerLine, lps)' - offset; 	% Converts data into proper shape for frame
				elseif state.acq.dualLaserMode==2	% alternate by line
					tempImageWhole = reshape(frameFinalData(:,inputChannelCounter), state.internal.samplesPerLine, 2*lps)' - offset; 	% Converts data into proper shape for frame
					tempImage = tempImageWhole(1:2:end, :);
					tempImage2 = tempImageWhole(2:2:end, :);
				end
				
				if state.acq.bidi
					tempImage(2:2:end,:)=fliplr(tempImage(2:2:end,:));
				end

				tempImage = add2d(tempImage(:, state.internal.startDataColumnInLine:state.internal.endDataColumnInLine), state.acq.binFactor);
				
				if state.acq.averaging & state.internal.frameCounter>1
					state.acq.acquiredData{channelCounter}(startLine:stopLine,:,position) = ...
						(...
						(state.internal.frameCounter - 1) ...	
						* state.acq.acquiredData{channelCounter}(startLine:stopLine,:,position...
						) ...
						+ tempImage ...
						)...
						/state.internal.frameCounter;
				else
					state.acq.acquiredData{channelCounter}(startLine:stopLine,:,position) = tempImage;
				end
				
				% Displays the current images on the screen as they are acquired.
				if state.acq.imagingChannel(channelCounter)
					set(state.internal.imagehandle(channelCounter), 'EraseMode', 'none', 'CData', ...
						state.acq.acquiredData{channelCounter}(startLine:stopLine,:,position), 'YData', [startLine stopLine]);
				end
				if state.internal.composite
					low=getfield(state.internal, ['lowPixelValue' num2str(channelCounter)]);
					hi=getfield(state.internal, ['highPixelValue' num2str(channelCounter)]);
					state.acq.compositeData(startLine:stopLine,:,state.internal.colorPlace(channelCounter))=...
						min(max(...
						(state.acq.acquiredData{channelCounter}(startLine:stopLine,:,position) - low) / max(hi-low, 1)...
						,0)...
						,1);
                end
				
				if state.acq.dualLaserMode==2	% alternate by line repeat for the other data
					if state.acq.bidi
						tempImage2(2:2:end,:)=fliplr(tempImage2(2:2:end,:));
					end
	
					tempImage2 = add2d(tempImage2(:, state.internal.startDataColumnInLine:state.internal.endDataColumnInLine), state.acq.binFactor);
					
					if state.acq.averaging & state.internal.frameCounter>1
						state.acq.acquiredData{10+channelCounter}(startLine:stopLine,:,position) = ...
							(...
							(state.internal.frameCounter - 1) ...	
							* state.acq.acquiredData{10+channelCounter}(startLine:stopLine,:,position...
							) ...
							+ tempImage2 ...
							)...
							/state.internal.frameCounter;
					else
						state.acq.acquiredData{10+channelCounter}(startLine:stopLine,:,position) = tempImage2;
					end
					
					% Displays the current images on the screen as they are acquired.
					if state.acq.imagingChannel(channelCounter)
						set(state.internal.imagehandle(10+channelCounter), 'EraseMode', 'none', 'CData', ...
							state.acq.acquiredData{10+channelCounter}(startLine:stopLine,:,position), 'YData', [startLine stopLine]);
					end
				end
				
			end
		end

		if state.internal.composite
			set(state.internal.compositeImagehandle, 'EraseMode', 'none', 'CData', ...
				state.acq.compositeData(startLine:stopLine,:,:), 'YData', [startLine stopLine]);
		end

		drawnow; 
		if state.internal.abortActionFunctions
			abortGrab;
			return
		end

		state.internal.stripeCounter = state.internal.stripeCounter + 1;
		%	disp([num2str(state.internal.stripeCounter) ' ' num2str(state.internal.frameCounter)]);
		
		if state.internal.stripeCounter == state.internal.numberOfStripes
			state.internal.stripeCounter = 0;
			if state.internal.frameCounter == state.acq.numberOfFrames 
				updateGUIByGlobal('state.internal.frameCounter');
				endAcquisition; 	% ResumeLoop, parkLaser, Close Shutter, appendData, reset counters,...
			end
		end
	catch
		if state.internal.abortActionFunctions
			abortGrab;
			return
		else
			setStatusString('Error in frame by stripes');
			disp('makeFrameByStripes: Error in action function');
			disp(lasterr);
		end
	end

		
	
	
