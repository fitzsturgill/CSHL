function preallocateMemory
	global state

	% This function Preallocates the appropriate memory for each acquisition mode.
	%
	% Written by: Thomas Pologruto and Bernardo Sabatini
	% Cold Spring Harbor Labs
	% February 1, 2001

    if state.acq.dualLaserMode==1   % lasers go simulataneously therefor one image window per laser
        channelList=1:state.init.maximumNumberOfInputChannels;
    else
        channelList=[1:state.init.maximumNumberOfInputChannels 11:10+state.init.maximumNumberOfInputChannels];
    end

	if ~iscell(state.acq.maxData)
		state.acq.maxData = cell(1,10+state.init.maximumNumberOfInputChannels);
	end
	state.acq.acquiredData = cell(1,10+state.init.maximumNumberOfInputChannels);
	state.acq.focusData = cell(1,10+state.init.maximumNumberOfInputChannels);

	
	for channelCounter = channelList
        inputChannelCounter=mod(channelCounter, 10);
		if getfield(state.acq, ['acquiringChannel' num2str(inputChannelCounter)])			% BSMOD 1/18/2 - removed eval for channelOn
			state.acq.focusData{channelCounter} = zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine);
			if state.acq.numberOfZSlices == 1 | state.internal.keepAllSlicesInMemory==0		% BSMOD 1/18/2 - added or statement
				% Continuous Time Series or only saving 1 slice in memory
				if state.acq.averaging == 0			% No averaging
					state.acq.acquiredData{channelCounter}=zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine, 1);
					state.acq.acquiredData{channelCounter}(:,:,state.acq.numberOfFrames)= ...
						zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine);
						
				elseif state.acq.averaging == 1 	% Averaging Time Series
					state.acq.acquiredData{channelCounter}= zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine, 1);
					
				end

			elseif state.acq.numberOfZSlices > 1 	% Discontinuous Z-Stack
				if state.acq.averaging == 0			% No averaging
					state.acq.acquiredData{channelCounter}=zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine, 1);
					state.acq.acquiredData{channelCounter}(:,:,state.acq.numberOfFrames*state.acq.numberOfZSlices) = ...
						zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine);
		
				elseif state.acq.averaging == 1 	% Averaging Z-Stack Series
					state.acq.acquiredData{channelCounter}= zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine, 1);
					state.acq.acquiredData{channelCounter}(:,:,state.acq.numberOfZSlices)= ...
						zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine);
			
				end
			end
			if getfield(state.acq, ['maxImage' num2str(inputChannelCounter)])
				if (size(state.acq.maxData{channelCounter}, 1) ~= state.acq.linesPerFrame) |  ...
						(size(state.acq.maxData{channelCounter}, 2) ~= state.acq.pixelsPerLine)
					state.acq.maxData{channelCounter} = zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine);
				end
			end
		else
			state.acq.acquiredData{channelCounter}=[];
		end			
	end

	state.acq.compositeData = (zeros(state.acq.linesPerFrame, state.acq.pixelsPerLine, 3)); 	% BSMOD 7/17/2

