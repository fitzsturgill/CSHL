function founderr=runLineScanAnalysis
	global state
	
	founderr=0;
	% calculate line scan
	
	state.analysis.deltax = state.acq.msPerLine*1000;
	
	state.analysis.acquiredData=cell(1,state.init.maximumNumberOfInputChannels);
	
	for channel=1:state.init.maximumNumberOfInputChannels
		if ~isempty(state.acq.acquiredData{channel}) & getfield(state.analysis, ['anaChannel' num2str(channel)])
			scanSize = size(state.acq.acquiredData{channel});
			if ndims(state.acq.acquiredData{channel})==3
				disp(['runFluorAnalysis: Channel #' num2str(channel) ' data has 3 dimensions.  Will reshape to 2 for line scan analysis']);
				tempData = permute(state.acq.acquiredData{channel}, [2 1 3]);
				state.analysis.acquiredData{channel} = tempData(:,:)';
			else
				state.analysis.acquiredData{channel} = state.acq.acquiredData{channel};
			end
		end
	end

	state.analysis.avgLineScan=calcLineScan(...
		state.analysis.acquiredData{state.analysis.avgLineScanChannel}, ...
		state.analysis.avgLineScanStart, state.analysis.avgLineScanEnd, 1);	% set smoothing to 1

	if ~iswave('avgLineScanWave')
		wave('avgLineScanWave', state.analysis.avgLineScan);
	else
		setWave('avgLineScanWave', 'data', state.analysis.avgLineScan);
	end
	
	% 	determine ROIs
	offset=getfield(state.acq, ['pmtOffsetChannel' num2str(state.analysis.avgLineScanChannel)])*state.acq.binFactor;
	offsetAmp=offset;
	if getfield(state.acq, ['pmtOffsetAutoSubtractChannel' num2str(state.analysis.avgLineScanChannel)])
		offset=0;
	end
	if ~iswave('offsetWave')
		wave('offsetWave', repmat(offset, 1, evalin('base', 'length(avgLineScanWave)')));
	else
		setWave('offsetWave', 'data', repmat(offset, 1, evalin('base', 'length(avgLineScanWave)')));
	end
	
	if  state.analysis.autosetROI
		[roix, roiy]=findPeaks(state.analysis.avgLineScan, state.analysis.numberOfROI, offset+offsetAmp/5, state.analysis.roiWidth);
			
		%findROI_LineScan(state.analysis.avgLineScan,  state.analysis.numberOfROI, state.analysis.roiWidth, ...
	%		offset, offset+offsetAmp/2, offset+offsetAmp/5);
		state.analysis.roiDefs = roix;
	else
		for counter=1:state.analysis.numberOfROI
			roix = state.analysis.roiDefs(:,1:2);
			roiy = state.analysis.avgLineScan(roix);
		end
	end
	
	if (size(roix, 1) < state.analysis.numberOfROI) | (size(roix, 2) ~= 2)
		founderr=1;
		beep;
		disp('*** Need to select proper number of 1D ROIs or autoset could not find sufficient ROIs ***');
		setStatusString('INCORRECT 1D ROIs');
		return
	end

	for counter=1:state.analysis.numberOfROI
		if counter<=size(roix, 1)
            if ~iswave(['roi_' num2str(counter) 'x'])
    			wave(['roi_' num2str(counter) 'x'], roix(counter, :)-1);
            else
    			setWave(['roi_' num2str(counter) 'x'], 'data', roix(counter, :)-1);
            end
                
            if ~iswave(['roi_' num2str(counter) 'y'])
			    wave(['roi_' num2str(counter) 'y'], roiy(counter, :));
            else
    			setWave(['roi_' num2str(counter) 'y'], 'data', roiy(counter, :));
            end
		else
            if iswave(['roi_' num2str(counter) 'x'])
    			setWave(['roi_' num2str(counter) 'x'], 'data', []);
            end                
            if iswave(['roi_' num2str(counter) 'y'])
    			setWave(['roi_' num2str(counter) 'y'], 'data', []);
            end                
		end			
	end
	
	counter=state.analysis.numberOfROI+1;
	while iswave(['roi_' num2str(counter) 'x'])
% 		setWave(['roi_' num2str(counter) 'x'], 'data', []);
%       if iswave(['roi_' num2str(counter) 'y'])
% 			setWave(['roi_' num2str(counter) 'y'], 'data', []);
%       end    
		evalin('base', ['roi_' num2str(counter) 'x.data=[];']);
        if iswave(['roi_' num2str(counter) 'y'])
			evalin('base', ['roi_' num2str(counter) 'y.data=[];']);
        end                
		counter=counter+1;
	end

	% 	get fluorescence profile
	state.analysis.roiFluorData=cell(1, state.init.maximumNumberOfInputChannels);
	
	numLines=[];
 	for counter=1:state.init.maximumNumberOfInputChannels
		if getfield(state.analysis, ['anaChannel' num2str(counter)])
			offset=0;
			if getfield(state.analysis, ['autosubOffset' num2str(counter)]) & ~getfield(state.acq, ['pmtOffsetAutoSubtractChannel' num2str(counter)])
				offset=getfield(state.acq, ['pmtOffsetChannel' num2str(counter)])*state.acq.binFactor;
			end

			state.analysis.roiFluorData{counter} = roiFluor(state.analysis.acquiredData{counter}, state.analysis.roiDefs, offset);
		end
	end

