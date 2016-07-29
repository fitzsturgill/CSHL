function founderr=runFrameAnalysis
	global state
	founderr=0;
	% calculate line scan
	
	state.analysis.deltax = state.acq.msPerLine*1000*state.acq.linesPerFrame;
	
	% check rois
	if (size(state.analysis.roiDefs2D, 1) < state.analysis.numberOfROI) | (size(state.analysis.roiDefs2D, 2) ~= 4)
		founderr=1;
		beep;
		disp('*** Need to select proper number of 2D ROIs ***');
		setStatusString('INCORRECT 2D ROIs');
		return
	end
	
	offset=getfield(state.acq, ['pmtOffsetChannel' num2str(state.analysis.avgLineScanChannel)])*state.acq.binFactor;
	offsetAmp=offset;
	if getfield(state.acq, ['pmtOffsetAutoSubtractChannel' num2str(state.analysis.avgLineScanChannel)])
		offset=0;
	end

	% 	get fluorescence profile
	state.analysis.roiFluorData=cell(1, state.init.maximumNumberOfInputChannels);
	state.analysis.roiBaseLines=cell(1, state.init.maximumNumberOfInputChannels);
	state.analysis.roiScanMeans=cell(1, state.init.maximumNumberOfInputChannels);
	
	numPoints=[];
 	for counter=1:state.init.maximumNumberOfInputChannels
		if getfield(state.analysis, ['anaChannel' num2str(counter)])
			offset=0;
			if getfield(state.analysis, ['autosubOffset' num2str(counter)]) & ~getfield(state.acq, ['pmtOffsetAutoSubtractChannel' num2str(counter)])
				offset=getfield(state.acq, ['pmtOffsetChannel' num2str(counter)])*state.acq.binFactor;
			end

			state.analysis.roiFluorData{counter} = roiFluor(state.acq.acquiredData{counter}, state.analysis.roiDefs2D(1:state.analysis.numberOfROI, :), offset);
			if isempty(numPoints)
				numPoints=size(state.analysis.roiFluorData{counter},2);
				baseFrameStart = max(1, floor(state.analysis.flBaseLineStart/state.analysis.deltax));
				if state.analysis.flBaseLineEnd < 1 
					baseFrameEnd = numPoints;
				else
					baseFrameEnd = min(numPoints, ceil(state.analysis.flBaseLineEnd/state.analysis.deltax));
				end
			end
				
			state.analysis.roiBaseLines{counter} = mean(...
				state.analysis.roiFluorData{counter}(:, baseFrameStart:baseFrameEnd), ...
				2);
			state.analysis.roiScanMeans{counter} = mean(...
				state.analysis.roiFluorData{counter}, ...
				2);
		end
	end

