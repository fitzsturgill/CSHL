function runFluorAnalysis(acqNumber)
	global state
	if nargin<1
		acqNumber=[];
	end
	
	if state.analysis.analysisMode==1
		return
	elseif state.analysis.analysisMode==2 | (state.analysis.analysisMode==4 & state.acq.lineScanCB)
		if state.acq.lineScanCB 
			founderr=runLineScanAnalysis;
		else
			disp('*** runFluorAnalysis : data is 2d.  Skipping line scan analysis');
			return
		end
	elseif state.analysis.analysisMode==3 | (state.analysis.analysisMode==4 & ~state.acq.lineScanCB)
		if ~state.acq.lineScanCB 
			founderr=runFrameAnalysis;
		else
			disp('*** runFluorAnalysis : data is 1d.  Skipping frame analysis');
			return
		end	
	end

	if founderr
		return
	end
	
	% calculate baselines and means
	state.analysis.roiBaseline=cell(1, state.init.maximumNumberOfInputChannels);
	state.analysis.roiMean=cell(1, state.init.maximumNumberOfInputChannels);
	state.analysis.roiPeak=cell(1, state.init.maximumNumberOfInputChannels);

	baseLineStartPoint = round(1 + state.analysis.flBaseLineStart/state.analysis.deltax);
	baseLineEndPoint = round(1 + state.analysis.flBaseLineEnd/state.analysis.deltax);
	
	for counter=1:state.init.maximumNumberOfInputChannels
		if getfield(state.analysis, ['anaChannel' num2str(counter)])
			if state.analysis.analysisMode==2 | (state.analysis.analysisMode==4 & state.acq.lineScanCB)
				nROI=min(size(state.analysis.roiDefs, 1), state.analysis.numberOfROI);
			elseif state.analysis.analysisMode==3 | (state.analysis.analysisMode==4 & ~state.acq.lineScanCB)
				nROI=min(size(state.analysis.roiDefs2D, 1), state.analysis.numberOfROI);
			end
			
			for roiCounter=1:nROI
				if state.analysis.flBaseLineEnd < 1 
					baseLineEndPoint = size(state.analysis.roiFluorData{counter},2);
				end
				
				state.analysis.roiBaseline{counter} = mean(...
					state.analysis.roiFluorData{counter}(:, baseLineStartPoint:baseLineEndPoint), ...
					2);
				state.analysis.roiMean{counter} = mean(...
					state.analysis.roiFluorData{counter}, ...
					2);
				
				% store baselines, means in waves
				if ~isempty(acqNumber)
					for roiCounter=1:length(state.analysis.roiBaseline{counter})
						if ~iswave([ROIScanName(counter, roiCounter) '_f0'])
							wave([ROIScanName(counter, roiCounter) '_f0'], 0, 'xscale', [1 1]);
						end
						eval(['global ' ROIScanName(counter, roiCounter) '_f0']);
						eval([ROIScanName(counter, roiCounter) '_f0(acqNumber)=state.analysis.roiBaseline{counter}(roiCounter);']);
						
						if ~iswave([ROIScanName(counter, roiCounter) '_favg'])
							wave([ROIScanName(counter, roiCounter) '_favg'], 0, 'xscale', [1 1]);
						end
						eval(['global ' ROIScanName(counter, roiCounter) '_favg']);
						eval([ROIScanName(counter, roiCounter) '_favg(acqNumber)=state.analysis.roiMean{counter}(roiCounter);']);
					end
				end
			end
		end
	end
		
	
	% 	do ratioing
 	for counter=1:state.init.maximumNumberOfInputChannels
		if getfield(state.analysis, ['anaChannel' num2str(counter)])
			ratioMode = getfield(state.analysis, ['ratioMode' num2str(counter)]);
			if ratioMode <=3
				% normalize to baseline data
				state.analysis.roiFluorData{counter} = state.analysis.roiFluorData{counter}./...
					repmat(state.analysis.roiBaseline{ratioMode}, 1, size(state.analysis.roiFluorData{counter},2));
			elseif ratioMode<=6
				% normalize to channel point by point
				state.analysis.roiFluorData{counter} = state.analysis.roiFluorData{counter}./state.analysis.roiFluorData{ratioMode-3};
			elseif ratioMode<=9
				% normalize to the channel mean
				state.analysis.roiFluorData{counter} = state.analysis.roiFluorData{counter}./...
					repmat(state.analysis.roiMean{ratioMode-6}, 1, size(state.analysis.roiFluorData{counter},2));
			elseif ratioMode~=10
				% error What mode is it?
				error('runFluorAnalysis: unknown ratioMode');
			end
		end
	end

	% recalculate baselines post ratioing
	state.analysis.roiBaselineR=cell(1, state.init.maximumNumberOfInputChannels);
	for counter=1:state.init.maximumNumberOfInputChannels
		if getfield(state.analysis, ['anaChannel' num2str(counter)])
			for roiCounter=nROI
				if getfield(state.analysis, ['ratioMode' num2str(counter)])<=9
					state.analysis.roiBaselineR{counter} = mean(...
						state.analysis.roiFluorData{counter}(:, baseLineStartPoint:baseLineEndPoint), ...
						2);
				else
					state.analysis.roiBaselineR{counter} = state.analysis.roiBaseline{counter};
				end
				
				% store baselines, means, and peaks in waves
				if ~isempty(acqNumber)
					for roiCounter=1:length(state.analysis.roiBaseline{counter})
						if ~iswave([ROIScanName(counter, roiCounter) '_rf0'])
							wave([ROIScanName(counter, roiCounter) '_rf0'], 0, 'xscale', [1 1]);
						end
						eval(['global ' ROIScanName(counter, roiCounter) '_rf0']);
						eval([ROIScanName(counter, roiCounter) '_rf0(acqNumber)=state.analysis.roiBaselineR{counter}(roiCounter);']);
					end
				end
			end
		end
	end
	
	% produce waves with results
 	for counter=1:state.init.maximumNumberOfInputChannels
		if getfield(state.analysis, ['anaChannel' num2str(counter)])
			for roiCounter=1:nROI
				if iswave(ROIScanName(counter, roiCounter))
					setWave(ROIScanName(counter, roiCounter), ...
						'data', state.analysis.roiFluorData{counter}(roiCounter, :), ...
						'xscale', [0 state.analysis.deltax]);
				else
					wave(ROIScanName(counter, roiCounter), state.analysis.roiFluorData{counter}(roiCounter, :), ...
						'xscale', [0 state.analysis.deltax]);
				end
			%	setWaveUserDataField(ROIScanName(counter, roiCounter), 'headerString', state.headerString);
				if state.analysis.analysisMode==2 | (state.analysis.analysisMode==4 & state.acq.lineScanCB)
					setWaveUserDataField(ROIScanName(counter, roiCounter), 'ROIDef', state.analysis.roiDefs(roiCounter,:));
				elseif state.analysis.analysisMode==3 | (state.analysis.analysisMode==4 & ~state.acq.lineScanCB)
					setWaveUserDataField(ROIScanName(counter, roiCounter), 'ROIDef', state.analysis.roiDefs2D(roiCounter,:));
				end

				if ~isempty(acqNumber)
					duplicateo(ROIScanName(counter, roiCounter), ROIScanName(counter, roiCounter, acqNumber));
				end
			end
		end
	end

	avginROIScans(acqNumber);
	if state.files.autoSave
		saveROIScans(acqNumber);
	end
	
	try
		if state.analysis.active
			runTraceAnalyzer(0);
		end
	catch
		disp(['runFluorAnalysis : ' lasterr]);
		disp('	when doing trace analysis');
	end
	
	if ~state.analysis.keepInMemory
	 	for counter=1:state.init.maximumNumberOfInputChannels
			for roiCounter=nROI
				kill(ROIScanName(counter, roiCounter, acqNumber));
			end
		end
	end
	
% 	